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
// $Id: SimpleExplicitIntegrator.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_SIMPLE_EXPLICIT_INTEGRATOR_H__
#define INCLUDED_SIMPLE_EXPLICIT_INTEGRATOR_H__

#include "Integrator.h"
#include "StructuredDataExchanger.h"
#include "IterationContext.h"

template <class Model, class ResEval>
class SimpleExplicitIntegrator : public Integrator
{
public:
    SimpleExplicitIntegrator(
        const Model& model,
        ResEval* resEval, // takes ownership
        Block& block, Structured<double>& U, const Structured<double>& U2, const Structured<double>& U3
        );
    virtual ~SimpleExplicitIntegrator();

    virtual const Structured<double>& ResidualVector() const { return mR; }

    virtual void PreIntegrateStart(int step, const IterationContext& iteration);
    virtual void PreIntegrateFinish(int step, const IterationContext& iteration);

    virtual void Integrate(
        int step, const IterationContext& iteration,
        const Structured<double>& U
        );

    virtual void PostIntegrateStart(int step, const IterationContext& iteration);
    virtual void PostIntegrateFinish(int step, const IterationContext& iteration);
    virtual const Structured<double>& RHS() const { return mR; }
    virtual Residual LatestResidual() const { return mResidual; }

protected:

private:
    Model mModel;
    ResEval* mResEval;

    Structured<double> mR;
    StructuredDataExchanger mSDX;

    Residual mResidual;
    const Structured<double>& mU2, mU3; // FIXME
};

#include "SimpleExplicitIntegrator.ipp"

#endif // INCLUDED_SIMPLE_EXPLICIT_INTEGRATOR_H__

