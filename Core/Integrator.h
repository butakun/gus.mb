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
// $Id: Integrator.h 255 2012-12-05 07:47:00Z kato $
#ifndef INCLUDED_INTEGRATOR_H__
#define INCLUDED_INTEGRATOR_H__

#include "Block.h"
#include "Structured.h"
#include "Residual.h"
#include <vector>

class IterationContext;

class Integrator
{
public:
    Integrator(int dof, Block& block, int numSteps);
    virtual ~Integrator();

    /// Returns the number of steps the integration needs, in order to complete one integration cycle.
    int IntegrationSteps() const { return mIntegrationSteps; }

    Structured<double>& DU() { return mDU; }
    const Structured<double>& DU() const { return mDU; }
    Structured<double>& DT() { return mDT; }
    const Structured<double>& DT() const { return mDT; }

    virtual const Structured<double>& ResidualVector() const = 0;

    virtual void PreIntegrateStart(int step, const IterationContext& iteration) = 0;
    virtual void PreIntegrateFinish(int step, const IterationContext& iteration) = 0;

    virtual void Integrate(
        int step, const IterationContext& iteration,
        const Structured<double>& U
        ) = 0;

    virtual void PostIntegrateStart(int step, const IterationContext& iteration) = 0;
    virtual void PostIntegrateFinish(int step, const IterationContext& iteration) = 0;

    virtual const Structured<double>& RHS() const = 0;
    virtual Residual LatestResidual() const = 0;

    Block& GetBlock() { return mBlock; }
    const Block& GetBlock() const { return mBlock; }

protected:

private:
    Block& mBlock;
    int mIntegrationSteps;
    Structured<double> mDU;
    Structured<double> mDT;
};

class Integrators : public std::vector<Integrator*>
{
public:
    Integrators() {}

    void PreIntegrateStart(int step, const IterationContext& iteration);
    void PreIntegrateFinish(int step, const IterationContext& iteration);

    void Integrate(
        int step, const IterationContext& iteration,
        const Structured<double>& U
        );

    void PostIntegrateStart(int step, const IterationContext& iteration);
    void PostIntegrateFinish(int step, const IterationContext& iteration);

protected:

private:
    Blocks mBlocks;
};

#endif // INCLUDED_INTEGRATOR_H__

