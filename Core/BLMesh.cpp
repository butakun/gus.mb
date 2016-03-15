/*
    gus.mb, an open source flow solver.
    Copyright (C) 2016 Hiromasa Kato <hiromasa at gmail.com>

    This file is part of gus.mb.

    gus.mb is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    gus.mb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
// $Id: BLMesh.cpp 123 2011-08-24 15:24:48Z kato $

#include "BLMesh.h"
#include "Vector3.h"
#include <cmath>
#include <cassert>

template <class T>
double Secant(T f, double x01, double x02, double ftol = 1.0e-5, int maxiter = 100)
{
    double x1, x2, x3 = 0.0;
    double f1, f2, f3;

    x1 = x01;
    x2 = x02;
    f1 = f(x1);
    f2 = f(x2);
    for (int i = 0; i < maxiter; ++i)
    {
        double dfdx = (f2 - f1) / (x2 - x1);
        // f = 0 = f2 + dfdx * (x3 - x2)
        double dx = -f2 / dfdx;
        for (int j = 0; j < 1000; ++j)
        {
            x3 = x2 + dx;
            if (x3 <= 1.0)
                dx *= 0.5;
        }
        f3 = f(x3);
        std::cout << i << ", " << x1 << ", " << f1 << ", " << x2 << ", " << f2 << ", " << x3 << ", " << f3 << std::endl;
        if (std::abs(f3) < ftol)
            break;
        x1 = x2; f1 = f2;
        x2 = x3; f2 = f3;
    }
    return x3;
}

double Roberts(double eta, double beta)
{
    double bp1 = beta + 1.0, bm1 = beta - 1.0;

    double A = std::pow(bp1 / bm1, 1.0 - eta);
    return (bp1 - bm1 * A) / (A + 1.0);
}

class RobertsFunc
{
public:
    RobertsFunc(double height_, double dywall_, int npts_) : height(height_), dywall(dywall_), npts(npts_) {}

    double operator () (double beta) const
    {
        double eta = 1.0 / double(npts - 1);
        return Roberts(eta, beta) - dywall / height;
    }

    double height, dywall;
    int npts;
};

class PowerFunc
{
public:
    PowerFunc(double Y_, double dy_, int npts_) : Y(Y_), dy(dy_), npts(npts_) {}

    double operator () (double a) const
    {
        double eta = 1.0 / double(npts - 1);
        return std::pow(eta, a) - dy / Y;
    }

    static double Func(double a, double eta)
    {
        return std::pow(eta, a);
    }

    double Y, dy;
    int npts;
};

void
BLMesh::GenerateMesh(Structured<double>& xyz, double LX, double LY, double LZ, int NLE, double XLE, double dywall) const
{
    assert(xyz.Data != NULL);

    bool bl = dywall > 0.0;

    IndexRange range = xyz.GetRange();

    double beta = 1.1;
    if (bl)
    {
        RobertsFunc rf(LY, dywall, range.End.J - range.Start.J + 1);
        beta = Secant<RobertsFunc>(rf, 1.1, 1.05);
        std::cout << "beta = " << beta << std::endl;
    }

    int ILE = range.Start.I + NLE - 1;
    double LXPLATE = LX - XLE;
    int NXPLATE = range.End.I - range.Start.I - NLE + 2;
    double alpha1 = 1.0, alpha2 = 1.0;
    std::cout << "NLE = " << NLE << ", NXPLATE = " << NXPLATE << ", ILE = " << ILE << std::endl;

    if (XLE > 0.0)
    {
        double dXLE = LXPLATE / double(NXPLATE - 1) * 0.1;
        std::cout << "dXLE = " << dXLE << std::endl;
        PowerFunc qf1(XLE, dXLE, NLE);
        alpha1 = Secant<PowerFunc>(qf1, 1.5, 2.0);
        std::cout << "alpha1 = " << alpha1 << ", dXLE = " << PowerFunc::Func(alpha1, 1.0 / double(NLE - 1)) * XLE << std::endl;
        PowerFunc qf2(LXPLATE, dXLE, NXPLATE);
        alpha2 = Secant<PowerFunc>(qf2, 1.5, 2.0);
        std::cout << "alpha2 = " << alpha2 << ", dXLE = " << PowerFunc::Func(alpha2, 1.0 / double(NXPLATE - 1)) * LXPLATE << std::endl;
    }

    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        double zeta = double(k - range.Start.K) / double(range.End.K - range.Start.K);
        double z = LZ * zeta;
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            bool debug = false && j == range.Start.J && k == range.Start.K;
            double eta = double(j - range.Start.J) / double(range.End.J - range.Start.J);
            double y;
            if (bl)
            {
                y = LY * Roberts(eta, beta);
            }
            else
            {
                y = LY * eta;
            }
            for (int i = range.Start.I; i <= ILE; ++i)
            {
                double xi = double(i - ILE) / double(range.Start.I - ILE);
                double x = XLE - XLE * PowerFunc::Func(alpha1, xi);
                if (debug) std::cout << i << ", " << xi << ", " << x << std::endl;
                xyz(i, j, k)[0] = x;
                xyz(i, j, k)[1] = y;
                xyz(i, j, k)[2] = z;
            }
            for (int i = ILE + 1; i <= range.End.I; ++i)
            {
                double xi = double(i - ILE) / double(range.End.I - ILE);
                double x = LXPLATE * PowerFunc::Func(alpha2, xi) + XLE;
                if (debug) std::cout << i << ", " << xi << ", " << x << std::endl;
                xyz(i, j, k)[0] = x;
                xyz(i, j, k)[1] = y;
                xyz(i, j, k)[2] = z;
            }
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                if (debug) std::cout << "i = " << i << ", X = " << xyz(i, j, k)[0] << std::endl;
            }
        }
    }
}

