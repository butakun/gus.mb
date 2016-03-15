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
// $Id: GMRESIntegrator.h 121 2011-08-18 16:08:26Z kato $
#ifndef INCLUDED_GMRES_INTEGRATOR_H__
#define INCLUDED_GMRES_INTEGRATOR_H__

#include "Integrator.h"
#include "LinearAlgebra.h"
#include "GMRES.h"

class BigVector
{
public:
    BigVector() : mU(NULL) {}
    ~BigVector()
    {
        if (mU != NULL)
        {
            delete[] mU->Data;
            delete mU;
        }
    }

    void Allocate(const Block& block)
    {
        if (mU != NULL)
        {
            delete[] mU->Data;
            delete mU;
        }
        mU = new Structured<double>(5, block.Dim());
    }

    Structured<double>& U() const { return *mU; }

    // X = a * X + b * Y, where X is myself and returns a reference to myself.
    BigVector& AXPBY(double a, double b, const BigVector& y)
    {
        // FIXME
        Structured<double>& X = *mU;
        Structured<double>& Y = y.U();
        IndexRange dim = mU->GetRange();
        for (int k = dim.Start.K; k <= dim.End.K; ++k)
            for (int j = dim.Start.J; j <= dim.End.J; ++j)
                for (int i = dim.Start.I; i <= dim.End.I; ++i)
                    for (int l = 0; l < mU->DOF(); ++l)
                        X(i, j, k)[l] = a * X(i, j, k)[l] + b * Y(i, j, k)[l];

        return *this;
    }

    double Norm() const
    {
        return std::sqrt(Dot(*this));
    }

    double Dot(const BigVector& v) const
    {
        const Structured<double>& V = v.U();
        double dot = 0.0;
        IndexRange dim = mU->GetRange();
        for (int k = dim.Start.K; k <= dim.End.K; ++k)
        {
            for (int j = dim.Start.J; j <= dim.End.J; ++j)
            {
                for (int i = dim.Start.I; i <= dim.End.I; ++i)
                {
                    for (int l = 0; l < mU->DOF(); ++l)
                    {
                        double u = (*mU)(i, j, k)[l];
                        double v = V(i, j, k)[l];
                        dot += u * v;
                    }
                }
            }
        }
        // FIXME: dot should be gathered from all the blocks
        return dot;
    }

    BigVector& operator = (double v)
    {
        mU->SetTo(v);
        return *this;
    }

    BigVector& operator * (double v)
    {
        mU->SetTo(v);
        return *this;
    }

protected:

private:
    Structured<double>* mU;
};

double norm(const BigVector& v) { return v.Norm(); }
double dot(const BigVector& v1, const BigVector& v2) { return v1.Dot(v2); }

class Preconditioner
{
public:
    Preconditioner(const Block& block) : mBlock(block) {}

    void Solve(BigVector& res, const BigVector& x) const
    {
        Structured<double>& Res = res.U();
        Structured<double>& X = x.U();
        IndexRange range = X.GetRange();
        for (int k = range.Start.K; k <= range.End.K; ++k)
        {
            for (int j = range.Start.J; j <= range.End.J; ++j)
            {
                for (int i = range.Start.I; i <= range.End.I; ++i)
                {
                    for (int l = 0; l < X.DOF(); ++l)
                    {
                        Res(i, j, k)[l] = X(i, j, k)[l];
                    }
                }
            }
        }
    }

protected:

private:
    const Block& mBlock;
};

class Ax
{
public:
    Ax(const Block& block, const Structured<double>& Rn);
    virtual ~Ax();

    BigVector& InitializeBigVector(BigVector& v) const;

    void Apply(BigVector& res, const BigVector& x) const;

protected:

private:
    const Structured<double>& mRn; // Residual at time n. (Remember we are solving dR/dU DU = -R)
};

template <class Model>
class GMRESIntegrator : public Integrator
{
public:
    GMRESIntegrator(const Model& model, Block& block);
    virtual ~GMRESIntegrator();

    virtual void PreIntegrateStart(int step);
    virtual void PreIntegrateFinish(int step);

    virtual void Integrate(
        int step,
        const Structured<double>& U
        );

    virtual void PostIntegrateStart(int step);
    virtual void PostIntegrateFinish(int step);

    virtual Residual LatestResidual() const;

    typedef GMRES<Ax, BigVector, Preconditioner, LinearAlgebra::Vector<double>, LinearAlgebra::Matrix<double>, double>
        GMRESSolver;

protected:

private:
    Model mModel;

    Structured<double> mRHS;

    Ax mAx;
    Preconditioner mPreconditioner;
    GMRESSolver mGMRES;
};

#include "GMRESIntegrator.ipp"

#endif // INCLUDED_GMRES_INTEGRATOR_H__
