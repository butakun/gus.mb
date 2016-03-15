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
// $Id: RectMesh.cpp 90 2011-01-28 14:19:59Z kato $

#include "RectMesh.h"
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

void
RectMesh::GenerateMesh(Structured<double>& xyz, double LX, double LY, double LZ, double dywall) const
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

    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        double zeta = double(k - range.Start.K) / double(range.End.K - range.Start.K);
        double z = LZ * zeta;
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
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
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                double xi = double(i - range.Start.I) / double(range.End.I - range.Start.I);
                double x = LX * xi;
                xyz(i, j, k)[0] = x;
                xyz(i, j, k)[1] = y;
                xyz(i, j, k)[2] = z;
            }
        }
    }
}

