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
// $Id: LUSGSIntegrator.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_LUSGS_INTEGRATOR_H__
#define INCLUDED_LUSGS_INTEGRATOR_H__

#include "Integrator.h"
#include "StructuredDataExchanger.h"
#include "Residual.h"
#include "IterationContext.h"

class ReducedFluxModel
{
public:
    ReducedFluxModel();
    double ScalarCoeff(double* Ui, double* Uj, double* Sij, double SijSign, double SijAbs) const;
    void JacobianDU(double* JacDU, double* U, double* dU, double* Sij, double SijSign, double SijAbs) const;

private:
    double Gamma;
};

template <class Model, class ResEval>
class LUSGSIntegrator : public Integrator
{
public:
    LUSGSIntegrator(
        const Model& model,
        ResEval* resEval, // FIXME: takes ownership
        Block& block,
        Structured<double>& U,
        const Structured<double>& U2,
        const Structured<double>& U3
        );
    virtual ~LUSGSIntegrator();

    virtual const Structured<double>& ResidualVector() const { return mRHS; }

    virtual void PreIntegrateStart(int step, const IterationContext& iteration);
    virtual void PreIntegrateFinish(int step, const IterationContext& iteration);
    virtual void Integrate(int step, const IterationContext& iteration, const Structured<double>& U);
    virtual void PostIntegrateStart(int step, const IterationContext& iteration);
    virtual void PostIntegrateFinish(int step, const IterationContext& iteration);
    virtual const Structured<double>& RHS() const { return mRHS; }
    virtual Residual LatestResidual() const { return mResidual; }

    const Structured<double>& DU2() const { return mDU2; }

protected:
    void LSweep(const Structured<double>& U, const IterationContext& iteration);
    void USweep(const Structured<double>& U);

    void USweep2(const Structured<double>& U, const IterationContext& iteration);
    void LSweep2(const Structured<double>& U);

    template <class Context> void Sweep1(const Structured<double>& U, const IterationContext& iteration);
    template <class Context> void Sweep2(const Structured<double>& U, const IterationContext& iteration);

private:
    //ReducedFluxModel mModel;
    Model mModel;
    ResEval* mResEval;

    Structured<double> mRHS; // RHS
    Structured<double> mRHS2; // FIXME: RHS (to keep the steady RHS for unsteady calc)
    Structured<double> mDU2; // L-U intermediate solution vector
    Structured<double> mD; // diagonal elements
    const Structured<double>& mU2; // FIXME: U at time level n-1
    const Structured<double>& mU3; // FIXME: U at time level n-2

    StructuredDataExchanger mSDXL;
    StructuredDataExchanger mSDXU;

    Residual mResidual;

    double mBDiag;

    bool mAlternateSweep;
};

#include "LUSGSIntegrator.ipp"

#endif // INCLUDED_LUSGS_INTEGRATOR_H__

