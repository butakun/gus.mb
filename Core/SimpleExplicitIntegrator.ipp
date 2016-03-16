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
// $Id: SimpleExplicitIntegrator.ipp 277 2013-06-04 01:58:51Z kato $

template <class Model, class ResEval>
SimpleExplicitIntegrator<Model, ResEval>::SimpleExplicitIntegrator(
    const Model& model, ResEval* resEval,
    Block& block, Structured<double>& U, const Structured<double>& U2, const Structured<double>& U3
    )
:   Integrator(mModel.DOF(), block, 1),
    mModel(model), mResEval(resEval),
    mR(model.DOF(), block.CellRange()),
    mSDX(block, U),
    mU2(U2), mU3(U3)
{
}

template <class Model, class ResEval>
SimpleExplicitIntegrator<Model, ResEval>::~SimpleExplicitIntegrator()
{
    delete[] mR.Data;
}

template <class Model, class ResEval>
void
SimpleExplicitIntegrator<Model, ResEval>::PreIntegrateStart(int step, const IterationContext& iteration)
{
}

template <class Model, class ResEval>
void
SimpleExplicitIntegrator<Model, ResEval>::PreIntegrateFinish(int step, const IterationContext& iteration)
{
}

template <class Model, class ResEval>
void
SimpleExplicitIntegrator<Model, ResEval>::Integrate(
    int step, const IterationContext& iteration,
    const Structured<double>& U
    )
{
    assert(step == 0);

    bool unsteady = iteration.IsUnsteady();
    double dTRealNonDim = iteration.DTNonDim();

    const Block& block = GetBlock();
    Structured<double>& dU = DU();
    Structured<double>& dT = DT();
    int dof = mModel.DOF();

    mResEval->EvaluateResidual(block, U, mR);
    mR.Multiply(-1.0);

    const Structured<double>& Vol = block.Vol();

    IndexRange cellRange = block.CellRange();
    for (int k = cellRange.Start.K; k <= cellRange.End.K; ++k)
    {
        for (int j = cellRange.Start.J; j <= cellRange.End.J; ++j)
        {
            for (int i = cellRange.Start.I; i <= cellRange.End.I; ++i)
            {
                double dt = *dT(i, j, k);
                double D = 1.0 / dt;
                double* rhs = mR(i, j, k);

                if (unsteady)
                {
                    double* un = mU2(i, j, k);
                    double* uk = U(i, j, k);
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] -= (uk[l] - un[l]) / dTRealNonDim;
                    }
                    D += 1.0 / dTRealNonDim;
                }

                for (int l = 0; l < dof; ++l)
                {
                    dU(i, j, k)[l] = 1.0 / D / *Vol(i, j, k) * rhs[l];
                }
            }
        }
    }

    mResidual = Residual(mR, cellRange, block.ID());
}

template <class Model, class ResEval>
void
SimpleExplicitIntegrator<Model, ResEval>::PostIntegrateStart(int step, const IterationContext& iteration)
{
    mModel.ApplyBCs(GetBlock());
    mSDX.Start();
}

template <class Model, class ResEval>
void
SimpleExplicitIntegrator<Model, ResEval>::PostIntegrateFinish(int step, const IterationContext& iteration)
{
    mSDX.Finish();
    GetBlock().FinalizeConnectivities(mSDX.Data());
}

