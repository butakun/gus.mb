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

#include "Vector3.h"
#include <cmath>
#include <iomanip>

template <class Model, class ResEval>
LUSGSIntegrator<Model, ResEval>::LUSGSIntegrator(
        const Model& model, ResEval* resEval,
        Block& block,
        Structured<double>& U,
        const Structured<double>& U2,
        const Structured<double>& U3
        )
:   Integrator(mModel.DOF(), block, 2),
    mModel(model), mResEval(resEval),
    mRHS(mModel.DOF(), block.CellRange()),
    mRHS2(mModel.DOF(), block.CellRange()),
    mDU2(mModel.DOF(), block.CellRangeWithGhosts()),
    mD(mModel.DOF(), block.CellRange()),
    mU2(U2), mU3(U3),
    mSDXL(block, &model, mDU2),
    mSDXU(block, &model, U),
    mBDiag(1.0),
    mAlternateSweep(false)
{
    mBDiag = 1.0;
}

template <class Model, class ResEval>
LUSGSIntegrator<Model, ResEval>::~LUSGSIntegrator()
{
    delete[] mRHS.Data;
    delete[] mRHS2.Data;
    delete[] mDU2.Data;
    delete[] mD.Data;
    delete mResEval;
}

template <class Model, class ResEval>
void
LUSGSIntegrator<Model, ResEval>::PreIntegrateStart(int step, const IterationContext& iteration)
{
}

template <class Model, class ResEval>
void
LUSGSIntegrator<Model, ResEval>::PreIntegrateFinish(int step, const IterationContext& iteration)
{
}

template <class Model, class ResEval>
void
LUSGSIntegrator<Model, ResEval>::Integrate(int step, const IterationContext& iteration, const Structured<double>& U)
{
    if (step == 0)
    {
        if (!mAlternateSweep)
        {
            LSweep(U, iteration);
        }
        else
        {
            USweep2(U, iteration);
        }
    }
    else
    {
        if (!mAlternateSweep)
        {
            USweep(U);
        }
        else
        {
            LSweep2(U);
        }
        mAlternateSweep = !mAlternateSweep;
    }
}

inline void DumpRHS(const Structured<double>& rhs, const IndexRange& range)
{
    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                std::cout << IndexIJK(i, j, k) << ": ";
                for (int l = 0; l < rhs.DOF(); ++l)
                {
                    std::cout << rhs(i, j, k)[l] << " ";
                }
                std::cout << std::endl;
            }
        }
    }
}

template <class Model, class ResEval>
void
LUSGSIntegrator<Model, ResEval>::LSweep(const Structured<double>& U, const IterationContext& iteration)
{
    bool unsteady = iteration.IsUnsteady();

    mResEval->EvaluateResidual(GetBlock(), U, mRHS);
    mRHS.Multiply(-1.0);

    double dTRealNonDim = iteration.DTNonDim();

    Structured<double>& dT = DT();

    const Block& block = GetBlock();
    IndexRange cr = block.CellRange();
    const Structured<double>& Sxi   = block.Sxi();
    const Structured<double>& Seta  = block.Seta();
    const Structured<double>& Szeta = block.Szeta();
    const Structured<double>& Vol = block.Vol();
    const Structured<double>& Radius = block.Radius();

    int dof = mModel.DOF();

    double* HdU = new double[dof];
    double* LdU = new double[dof];

    mDU2.SetTo(0.0);

    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* D = mD(i, j, k);

                bool xiM, etaM, zetaM;
                xiM   = i > cr.Start.I;
                etaM  = j > cr.Start.J;
                zetaM = k > cr.Start.K;

                double SxiPAbs, SxiMAbs, SetaPAbs, SetaMAbs, SzetaPAbs, SzetaMAbs;
                SxiPAbs = Vector3(Sxi(i, j, k)).Mag();
                SxiMAbs = Vector3(Sxi(i - 1, j, k)).Mag();
                SetaPAbs = Vector3(Seta(i, j, k)).Mag();
                SetaMAbs = Vector3(Seta(i, j - 1, k)).Mag();
                SzetaPAbs = Vector3(Szeta(i, j, k)).Mag();
                SzetaMAbs = Vector3(Szeta(i, j, k - 1)).Mag();

                double dt = *dT(i, j, k);
                double vol = *Vol(i, j, k);
#define SCALAR_DIAG 0
#if SCALAR_DIAG
                if (dt > 0.0)
                    *mD(i, j, k) = vol / dt;
                else
                    *mD(i, j, k) = 0.0;
#else
                double dcoef = dt > 0.0 ? 1.0 : 0;
                for (int l = 0; l < dof; ++l)
                    D[l] = dcoef * vol / dt;
#endif

                IndexIJK i0, iXiP, iXiM, iEtaP, iEtaM, iZetaP, iZetaM;
                i0     = IndexIJK(i, j, k);
                iXiP   = IndexIJK(i + 1, j, k); iXiM   = IndexIJK(i - 1, j, k);
                iEtaP  = IndexIJK(i, j + 1, k); iEtaM  = IndexIJK(i, j - 1, k);
                iZetaP = IndexIJK(i, j, k + 1); iZetaM = IndexIJK(i, j, k - 1);

                double lamXiP, lamXiM, lamEtaP, lamEtaM, lamZetaP, lamZetaM;
                lamXiP   = mModel.ScalarCoeff(block, U, i0, iXiP,     Sxi(i,     j, k),  1.0, SxiPAbs, Radius);
                lamXiM   = mModel.ScalarCoeff(block, U, i0, iXiM,     Sxi(i - 1, j, k), -1.0, SxiMAbs, Radius);
                lamEtaP  = mModel.ScalarCoeff(block, U, i0, iEtaP,   Seta(i,     j, k),  1.0, SetaPAbs, Radius);
                lamEtaM  = mModel.ScalarCoeff(block, U, i0, iEtaM,   Seta(i, j - 1, k), -1.0, SetaMAbs, Radius);
                lamZetaP = mModel.ScalarCoeff(block, U, i0, iZetaP, Szeta(i, j, k    ),  1.0, SzetaPAbs, Radius);
                lamZetaM = mModel.ScalarCoeff(block, U, i0, iZetaM, Szeta(i, j, k - 1), -1.0, SzetaMAbs, Radius);

#if SCALAR_DIAG
                *mD(i, j, k) += 0.5 * (lamXiP + lamXiM + lamEtaP + lamEtaM + lamZetaP + lamZetaM);
                *mD(i, j, k) += mModel.DiagonalAddition(block, U, i0);
#else
                double d = 0.5 * (lamXiP + lamXiM + lamEtaP + lamEtaM + lamZetaP + lamZetaM);
                for (int l = 0; l < dof; ++l)
                    D[l] += d;
                mModel.DiagonalAddition(D, block, U, i0);
#endif

                for (int l = 0; l < dof; ++l)
                {
                    LdU[l] = 0.0;
                }
                if (xiM)
                {
                    mModel.JacobianDU(HdU, block, iXiM, U(iXiM), mDU2(iXiM), Sxi(i - 1, j, k), -1.0, SxiMAbs, Radius(iXiM));
                    for (int l = 0; l < dof; ++l)
                    {
                        LdU[l] += 0.5 * (HdU[l] - lamXiM * mDU2(iXiM)[l]);
                    }
                }
                if (etaM)
                {
                    mModel.JacobianDU(HdU, block, iEtaM, U(iEtaM), mDU2(iEtaM), Seta(i, j - 1, k), -1.0, SetaMAbs, Radius(iEtaM));
                    for (int l = 0; l < dof; ++l)
                    {
                        LdU[l] += 0.5 * (HdU[l] - lamEtaM * mDU2(iEtaM)[l]);
                    }
                }
                if (zetaM)
                {
                    mModel.JacobianDU(HdU, block, iZetaM, U(iZetaM), mDU2(iZetaM), Szeta(i, j, k - 1), -1.0, SzetaMAbs, Radius(iZetaM));
                    for (int l = 0; l < dof; ++l)
                    {
                        LdU[l] += 0.5 * (HdU[l] - lamZetaM * mDU2(iZetaM)[l]);
                    }
                }

                if (!xiM || !etaM || !zetaM)
                {
#if SCALAR_DIAG
                    *mD(i, j, k) *= mBDiag;
#else
                    for (int l = 0; l < dof; ++l)
                        D[l] *= mBDiag;
#endif
                }

                if (unsteady)
                {
                    double* uk = U(i, j, k);
                    double* un1 = mU2(i, j, k); // n
                    double* un2 = mU3(i, j, k); // n - 1
                    double* rhs = mRHS(i, j, k);

#define FIRST_ORDER 0
#if FIRST_ORDER
#if SCALAR_DIAG
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] += -vol * (uk[l] - un1[l]) / dTRealNonDim;
                    }
                    *mD(i, j, k) += vol / dTRealNonDim;
#else
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] += -vol * (uk[l] - un1[l]) / dTRealNonDim;
                        D[l] += vol / dTRealNonDim;
                    }
#endif
#else
#if SCALAR_DIAG
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] += -vol * 0.5 * (3.0 * uk[l] - 4.0 * un1[l] + un2[l]) / dTRealNonDim;
                    }
                    *mD(i, j, k) += 1.5 * vol / dTRealNonDim;
#else
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] += -vol * 0.5 * (3.0 * uk[l] - 4.0 * un1[l] + un2[l]) / dTRealNonDim;
                        D[l] += 1.5 * vol / dTRealNonDim;
                    }
#endif
#endif
                }

                for (int l = 0; l < dof; ++l)
                {
#if SCALAR_DIAG
                    mDU2(i, j, k)[l] = (mRHS(i, j, k)[l] - LdU[l]) / (*mD(i, j, k));
#else
                    mDU2(i, j, k)[l] = (mRHS(i, j, k)[l] - LdU[l]) / D[l];
#endif
                }
            }
        }
    }

    mResidual = Residual(mRHS, block.CellRange(), block.ID());

    delete[] HdU;
    delete[] LdU;
}

template <class Model, class ResEval>
void
LUSGSIntegrator<Model, ResEval>::USweep(const Structured<double>& U)
{
    Structured<double>& dU = DU();

    const Block& block = GetBlock();
    IndexRange cr = block.CellRange();
    const Structured<double>& Sxi   = block.Sxi();
    const Structured<double>& Seta  = block.Seta();
    const Structured<double>& Szeta = block.Szeta();
    const Structured<double>& Radius = block.Radius();

    int dof = mModel.DOF();
    double* HdU = new double[dof];
    double* UdU = new double[dof];

    for (int k = cr.End.K; k >= cr.Start.K; --k)
    {
        for (int j = cr.End.J; j >= cr.Start.J; --j)
        {
            for (int i = cr.End.I; i >= cr.Start.I; --i)
            {
                double* D = mD(i, j, k);

                for (int l = 0; l < dof; ++l)
                {
#if SCALAR_DIAG
                    mDU2(i, j, k)[l] *= (*mD(i, j, k));
#else
                    mDU2(i, j, k)[l] *= D[l];
#endif
                }

                bool xiP, etaP, zetaP;
                xiP   = i < cr.End.I;
                etaP  = j < cr.End.J;
                zetaP = k < cr.End.K;

                double SxiPAbs, SetaPAbs, SzetaPAbs;
                SxiPAbs = Vector3(Sxi(i, j, k)).Mag();
                SetaPAbs = Vector3(Seta(i, j, k)).Mag();
                SzetaPAbs = Vector3(Szeta(i, j, k)).Mag();

                IndexIJK i0, iXiP, iEtaP, iZetaP;
                i0 = IndexIJK(i, j, k);
                iXiP   = IndexIJK(i + 1, j, k);
                iEtaP  = IndexIJK(i, j + 1, k);
                iZetaP = IndexIJK(i, j, k + 1);

                double lamXiP, lamEtaP, lamZetaP;
                lamXiP   = mModel.ScalarCoeff(block, U, i0, iXiP,     Sxi(i, j, k), 1.0, SxiPAbs, Radius);
                lamEtaP  = mModel.ScalarCoeff(block, U, i0, iEtaP,   Seta(i, j, k), 1.0, SetaPAbs, Radius);
                lamZetaP = mModel.ScalarCoeff(block, U, i0, iZetaP, Szeta(i, j, k), 1.0, SzetaPAbs, Radius);

                for (int l = 0; l < dof; ++l)
                {
                    UdU[l] = 0.0;
                }
                if (xiP)
                {
                    mModel.JacobianDU(HdU, block, iXiP, U(iXiP), mDU2(iXiP), Sxi(i, j, k), 1.0, SxiPAbs, Radius(iXiP));
                    for (int l = 0; l < dof; ++l)
                    {
                        UdU[l] += 0.5 * (HdU[l] - lamXiP * mDU2(iXiP)[l]);
                    }
                }
                if (etaP)
                {
                    mModel.JacobianDU(HdU, block, iEtaP, U(iEtaP), mDU2(iEtaP), Seta(i, j, k), 1.0, SetaPAbs, Radius(iEtaP));
                    for (int l = 0; l < dof; ++l)
                    {
                        UdU[l] += 0.5 * (HdU[l] - lamEtaP * mDU2(iEtaP)[l]);
                    }
                }
                if (zetaP)
                {
                    mModel.JacobianDU(HdU, block, iZetaP, U(iZetaP), mDU2(iZetaP), Szeta(i, j, k), 1.0, SzetaPAbs, Radius(iZetaP));
                    for (int l = 0; l < dof; ++l)
                    {
                        UdU[l] += 0.5 * (HdU[l] - lamZetaP * mDU2(iZetaP)[l]);
                    }
                }

                for (int l = 0; l < dof; ++l)
                {
#if SCALAR_DIAG
                    dU(i, j, k)[l] = (mDU2(i, j, k)[l] - UdU[l]) / (*mD(i, j, k));
#else
                    dU(i, j, k)[l] = (mDU2(i, j, k)[l] - UdU[l]) / D[l];
#endif
                }
            }
        }
    }

    delete[] HdU;
    delete[] UdU;
}

template <class Model, class ResEval>
void
LUSGSIntegrator<Model, ResEval>::USweep2(const Structured<double>& U, const IterationContext& iteration)
{
    bool unsteady = iteration.IsUnsteady();

    mResEval->EvaluateResidual(GetBlock(), U, mRHS);
    mRHS.Multiply(-1.0);

    double dTRealNonDim = iteration.DTNonDim();

    Structured<double>& dT = DT();

    const Block& block = GetBlock();
    IndexRange cr = block.CellRange();
    const Structured<double>& Sxi   = block.Sxi();
    const Structured<double>& Seta  = block.Seta();
    const Structured<double>& Szeta = block.Szeta();
    const Structured<double>& Vol = block.Vol();
    const Structured<double>& Radius = block.Radius();

    int dof = mModel.DOF();

    double* HdU = new double[dof];
    double* UdU = new double[dof];

    mDU2.SetTo(0.0);

    for (int k = cr.End.K; k >= cr.Start.K; --k)
    {
        for (int j = cr.End.J; j >= cr.Start.J; --j)
        {
            for (int i = cr.End.I; i >= cr.Start.I; --i)
            {
                double* D = mD(i, j, k);

                bool xiP, etaP, zetaP;
                xiP   = i < cr.End.I;
                etaP  = j < cr.End.J;
                zetaP = k < cr.End.K;

                double SxiPAbs, SxiMAbs, SetaPAbs, SetaMAbs, SzetaPAbs, SzetaMAbs;
                SxiPAbs = Vector3(Sxi(i, j, k)).Mag();
                SxiMAbs = Vector3(Sxi(i - 1, j, k)).Mag();
                SetaPAbs = Vector3(Seta(i, j, k)).Mag();
                SetaMAbs = Vector3(Seta(i, j - 1, k)).Mag();
                SzetaPAbs = Vector3(Szeta(i, j, k)).Mag();
                SzetaMAbs = Vector3(Szeta(i, j, k - 1)).Mag();

                double dt = *dT(i, j, k);
                double vol = *Vol(i, j, k);
#define SCALAR_DIAG 0
#if SCALAR_DIAG
                if (dt > 0.0)
                    *mD(i, j, k) = vol / dt;
                else
                    *mD(i, j, k) = 0.0;
#else
                double dcoef = dt > 0.0 ? 1.0 : 0;
                for (int l = 0; l < dof; ++l)
                    D[l] = dcoef * vol / dt;
#endif

                IndexIJK i0, iXiP, iXiM, iEtaP, iEtaM, iZetaP, iZetaM;
                i0     = IndexIJK(i, j, k);
                iXiP   = IndexIJK(i + 1, j, k); iXiM   = IndexIJK(i - 1, j, k);
                iEtaP  = IndexIJK(i, j + 1, k); iEtaM  = IndexIJK(i, j - 1, k);
                iZetaP = IndexIJK(i, j, k + 1); iZetaM = IndexIJK(i, j, k - 1);

                double lamXiP, lamXiM, lamEtaP, lamEtaM, lamZetaP, lamZetaM;
                lamXiP   = mModel.ScalarCoeff(block, U, i0, iXiP,     Sxi(i,     j, k),  1.0, SxiPAbs, Radius);
                lamXiM   = mModel.ScalarCoeff(block, U, i0, iXiM,     Sxi(i - 1, j, k), -1.0, SxiMAbs, Radius);
                lamEtaP  = mModel.ScalarCoeff(block, U, i0, iEtaP,   Seta(i,     j, k),  1.0, SetaPAbs, Radius);
                lamEtaM  = mModel.ScalarCoeff(block, U, i0, iEtaM,   Seta(i, j - 1, k), -1.0, SetaMAbs, Radius);
                lamZetaP = mModel.ScalarCoeff(block, U, i0, iZetaP, Szeta(i, j, k    ),  1.0, SzetaPAbs, Radius);
                lamZetaM = mModel.ScalarCoeff(block, U, i0, iZetaM, Szeta(i, j, k - 1), -1.0, SzetaMAbs, Radius);

#if SCALAR_DIAG
                *mD(i, j, k) += 0.5 * (lamXiP + lamXiM + lamEtaP + lamEtaM + lamZetaP + lamZetaM);
                *mD(i, j, k) += mModel.DiagonalAddition(block, U, i0);
#else
                double d = 0.5 * (lamXiP + lamXiM + lamEtaP + lamEtaM + lamZetaP + lamZetaM);
                for (int l = 0; l < dof; ++l)
                    D[l] += d;
                mModel.DiagonalAddition(D, block, U, i0);
#endif

                for (int l = 0; l < dof; ++l)
                {
                    UdU[l] = 0.0;
                }
                if (xiP)
                {
                    mModel.JacobianDU(HdU, block, iXiP, U(iXiP), mDU2(iXiP), Sxi(i, j, k), 1.0, SxiPAbs, Radius(iXiP));
                    for (int l = 0; l < dof; ++l)
                    {
                        UdU[l] += 0.5 * (HdU[l] - lamXiP * mDU2(iXiP)[l]);
                    }
                }
                if (etaP)
                {
                    mModel.JacobianDU(HdU, block, iEtaP, U(iEtaP), mDU2(iEtaP), Seta(i, j, k), 1.0, SetaPAbs, Radius(iEtaP));
                    for (int l = 0; l < dof; ++l)
                    {
                        UdU[l] += 0.5 * (HdU[l] - lamEtaP * mDU2(iEtaP)[l]);
                    }
                }
                if (zetaP)
                {
                    mModel.JacobianDU(HdU, block, iZetaP, U(iZetaP), mDU2(iZetaP), Szeta(i, j, k), 1.0, SzetaPAbs, Radius(iZetaP));
                    for (int l = 0; l < dof; ++l)
                    {
                        UdU[l] += 0.5 * (HdU[l] - lamZetaP * mDU2(iZetaP)[l]);
                    }
                }

                if (!xiP || !etaP || !zetaP)
                {
#if SCALAR_DIAG
                    *mD(i, j, k) *= mBDiag;
#else
                    for (int l = 0; l < dof; ++l)
                        D[l] *= mBDiag;
#endif
                }

                if (unsteady)
                {
                    double* uk = U(i, j, k);
                    double* un1 = mU2(i, j, k); // n
                    double* un2 = mU3(i, j, k); // n - 1
                    double* rhs = mRHS(i, j, k);

#if FIRST_ORDER
#if SCALAR_DIAG
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] += -vol * (uk[l] - un1[l]) / dTRealNonDim;
                    }
                    *mD(i, j, k) += vol / dTRealNonDim;
#else
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] += -vol * (uk[l] - un1[l]) / dTRealNonDim;
                        D[l] += vol / dTRealNonDim;
                    }
#endif
#else
#if SCALAR_DIAG
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] += -vol * 0.5 * (3.0 * uk[l] - 4.0 * un1[l] + un2[l]) / dTRealNonDim;
                    }
                    *mD(i, j, k) += 1.5 * vol / dTRealNonDim;
#else
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] += -vol * 0.5 * (3.0 * uk[l] - 4.0 * un1[l] + un2[l]) / dTRealNonDim;
                        D[l] += 1.5 * vol / dTRealNonDim;
                    }
#endif
#endif
                }

#if SCALAR_DIAG
                for (int l = 0; l < dof; ++l)
                {
                    mDU2(i, j, k)[l] = (mRHS(i, j, k)[l] - UdU[l]) / (*mD(i, j, k));
                }
#else
                for (int l = 0; l < dof; ++l)
                {
                    mDU2(i, j, k)[l] = (mRHS(i, j, k)[l] - UdU[l]) / D[l];
                }
#endif
            }
        }
    }

    mResidual = Residual(mRHS, block.CellRange(), block.ID());

    delete[] HdU;
    delete[] UdU;
}

template <class Model, class ResEval>
void
LUSGSIntegrator<Model, ResEval>::LSweep2(const Structured<double>& U)
{
    Structured<double>& dU = DU();

    const Block& block = GetBlock();
    IndexRange cr = block.CellRange();
    const Structured<double>& Sxi   = block.Sxi();
    const Structured<double>& Seta  = block.Seta();
    const Structured<double>& Szeta = block.Szeta();
    const Structured<double>& Radius = block.Radius();

    int dof = mModel.DOF();

    double* HdU = new double[dof];
    double* LdU = new double[dof];

    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* D = mD(i, j, k);

                for (int l = 0; l < dof; ++l)
                {
#if SCALAR_DIAG
                    mDU2(i, j, k)[l] *= (*mD(i, j, k));
#else
                    mDU2(i, j, k)[l] *= D[l];
#endif
                }

                bool xiM, etaM, zetaM;
                xiM   = i > cr.Start.I;
                etaM  = j > cr.Start.J;
                zetaM = k > cr.Start.K;

                double SxiMAbs, SetaMAbs, SzetaMAbs;
                SxiMAbs = Vector3(Sxi(i - 1, j, k)).Mag();
                SetaMAbs = Vector3(Seta(i, j - 1, k)).Mag();
                SzetaMAbs = Vector3(Szeta(i, j, k - 1)).Mag();

                IndexIJK i0, iXiP, iXiM, iEtaP, iEtaM, iZetaP, iZetaM;
                i0     = IndexIJK(i, j, k);
                iXiP   = IndexIJK(i + 1, j, k); iXiM   = IndexIJK(i - 1, j, k);
                iEtaP  = IndexIJK(i, j + 1, k); iEtaM  = IndexIJK(i, j - 1, k);
                iZetaP = IndexIJK(i, j, k + 1); iZetaM = IndexIJK(i, j, k - 1);

                double lamXiM, lamEtaM, lamZetaM;
                lamXiM   = mModel.ScalarCoeff(block, U, i0, iXiM,     Sxi(i - 1, j, k), -1.0, SxiMAbs, Radius);
                lamEtaM  = mModel.ScalarCoeff(block, U, i0, iEtaM,   Seta(i, j - 1, k), -1.0, SetaMAbs, Radius);
                lamZetaM = mModel.ScalarCoeff(block, U, i0, iZetaM, Szeta(i, j, k - 1), -1.0, SzetaMAbs, Radius);

                for (int l = 0; l < dof; ++l)
                {
                    LdU[l] = 0.0;
                }
                if (xiM)
                {
                    mModel.JacobianDU(HdU, block, iXiM, U(iXiM), mDU2(iXiM), Sxi(i - 1, j, k), -1.0, SxiMAbs, Radius(iXiM));
                    for (int l = 0; l < dof; ++l)
                    {
                        LdU[l] += 0.5 * (HdU[l] - lamXiM * mDU2(iXiM)[l]);
                    }
                }
                if (etaM)
                {
                    mModel.JacobianDU(HdU, block, iEtaM, U(iEtaM), mDU2(iEtaM), Seta(i, j - 1, k), -1.0, SetaMAbs, Radius(iEtaM));
                    for (int l = 0; l < dof; ++l)
                    {
                        LdU[l] += 0.5 * (HdU[l] - lamEtaM * mDU2(iEtaM)[l]);
                    }
                }
                if (zetaM)
                {
                    mModel.JacobianDU(HdU, block, iZetaM, U(iZetaM), mDU2(iZetaM), Szeta(i, j, k - 1), -1.0, SzetaMAbs, Radius(iZetaM));
                    for (int l = 0; l < dof; ++l)
                    {
                        LdU[l] += 0.5 * (HdU[l] - lamZetaM * mDU2(iZetaM)[l]);
                    }
                }

                for (int l = 0; l < dof; ++l)
                {
#if SCALAR_DIAG
                    dU(i, j, k)[l] = (mDU2(i, j, k)[l] - LdU[l]) / (*mD(i, j, k));
#else
                    dU(i, j, k)[l] = (mDU2(i, j, k)[l] - LdU[l]) / D[l];
#endif
                }
            }
        }
    }

    delete[] HdU;
    delete[] LdU;
}

template <class Model, class ResEval>
void
LUSGSIntegrator<Model, ResEval>::PostIntegrateStart(int step, const IterationContext& iteration)
{
    if (step == 0)
    {
        mSDXL.Start();
    }
    else if (step == 1)
    {
        //mModel.ApplyBCs(GetBlock());
        mSDXU.Start();
    }
    else
    {
        assert(false);
    }
}

template <class Model, class ResEval>
void
LUSGSIntegrator<Model, ResEval>::PostIntegrateFinish(int step, const IterationContext& iteration)
{
    if (step == 0)
    {
        mSDXL.Finish();
        GetBlock().FinalizeConnectivities(mDU2); // FIXME: shouldn't this mDU2 be accessed as mSDXL.Data()?
    }
    else if (step == 1)
    {
        mSDXU.Finish();
        GetBlock().FinalizeConnectivities(mSDXU.Data());
    }
    else
    {
        assert(false);
    }
}

#if 0
// Templatized sweeps

class SweepContext
{
public:
    IndexIJK Start1, End1, Increment1;
    IndexIJK Start2, End2, Increment2;

    IndexIJK IndexXiP(const IndexIJK& ijk) const;
    IndexIJK IndexXiM(const IndexIJK& ijk) const;
    IndexIJK IndexEtaP(const IndexIJK& ijk) const;
    IndexIJK IndexEtaM(const IndexIJK& ijk) const;
    IndexIJK IndexZetaP(const IndexIJK& ijk) const;
    IndexIJK IndexZetaM(const IndexIJK& ijk) const;
};

template <class Model, class ResEval>
template <class Context>
void
LUSGSIntegrator<Model, ResEval>::Sweep1<Context>(const Structured<double>& U, const IterationContext& iteration)
{
    bool unsteady = iteration.IsUnsteady();

    mModel.GetResidualEvaluator().EvaluateResidual(GetBlock(), U, mRHS);
    mRHS.Multiply(-1.0);

    double dTReal = iteration.DT();

    Structured<double>& dT = DT();

    const Block& block = GetBlock();
    IndexRange cr = block.CellRange();
    const Structured<double>& Sxi   = block.Sxi();
    const Structured<double>& Seta  = block.Seta();
    const Structured<double>& Szeta = block.Szeta();
    const Structured<double>& Vol = block.Vol();
    const Structured<double>& Radius = block.Radius();

    int dof = mModel.DOF();

    double* HdU = new double[dof];
    double* LdU = new double[dof];

    mDU2.SetTo(0.0);

    Context sweep(cr);

    for (int k = sweep.Start1.K; k <= sweep.End1.K; i += sweep.Increment1.I)
    {
        for (int j = sweep.Start1.J; j <= sweep.End1.J; j += sweep.Increment1.J)
        {
            for (int i = sweep.Start1.I; i <= sweep.End1.I; k += sweep.Increment1.K)
            {
                // FIXME
                bool xiM, etaM, zetaM;
                xiM   = i > cr.Start.I;
                etaM  = j > cr.Start.J;
                zetaM = k > cr.Start.K;

                double SxiPAbs, SxiMAbs, SetaPAbs, SetaMAbs, SzetaPAbs, SzetaMAbs;
                SxiPAbs = Vector3(Sxi(i, j, k)).Mag();
                SxiMAbs = Vector3(Sxi(i - 1, j, k)).Mag();
                SetaPAbs = Vector3(Seta(i, j, k)).Mag();
                SetaMAbs = Vector3(Seta(i, j - 1, k)).Mag();
                SzetaPAbs = Vector3(Szeta(i, j, k)).Mag();
                SzetaMAbs = Vector3(Szeta(i, j, k - 1)).Mag();

                double dt = *dT(i, j, k);
                double vol = *Vol(i, j, k);
                if (dt > 0.0)
                    *mD(i, j, k) = vol / dt;
                else
                    *mD(i, j, k) = 0.0;

                IndexIJK i0, iXiP, iXiM, iEtaP, iEtaM, iZetaP, iZetaM;
                i0     = IndexIJK(i, j, k);
                iXiP   = IndexIJK(i + 1, j, k); iXiM   = IndexIJK(i - 1, j, k);
                iEtaP  = IndexIJK(i, j + 1, k); iEtaM  = IndexIJK(i, j - 1, k);
                iZetaP = IndexIJK(i, j, k + 1); iZetaM = IndexIJK(i, j, k - 1);

                double lamXiP, lamXiM, lamEtaP, lamEtaM, lamZetaP, lamZetaM;
                lamXiP   = mModel.ScalarCoeff(block, U, i0, iXiP,     Sxi(i,     j, k),  1.0, SxiPAbs, Radius);
                lamXiM   = mModel.ScalarCoeff(block, U, i0, iXiM,     Sxi(i - 1, j, k), -1.0, SxiMAbs, Radius);
                lamEtaP  = mModel.ScalarCoeff(block, U, i0, iEtaP,   Seta(i,     j, k),  1.0, SetaPAbs, Radius);
                lamEtaM  = mModel.ScalarCoeff(block, U, i0, iEtaM,   Seta(i, j - 1, k), -1.0, SetaMAbs, Radius);
                lamZetaP = mModel.ScalarCoeff(block, U, i0, iZetaP, Szeta(i, j, k    ),  1.0, SzetaPAbs, Radius);
                lamZetaM = mModel.ScalarCoeff(block, U, i0, iZetaM, Szeta(i, j, k - 1), -1.0, SzetaMAbs, Radius);

                *mD(i, j, k) += 0.5 * (lamXiP + lamXiM + lamEtaP + lamEtaM + lamZetaP + lamZetaM);
                *mD(i, j, k) += mModel.DiagonalAddition(block, U, i0);

                for (int l = 0; l < dof; ++l)
                {
                    LdU[l] = 0.0;
                }
                if (xiM)
                {
                    mModel.JacobianDU(HdU, block, iXiM, U(iXiM), mDU2(iXiM), Sxi(i - 1, j, k), -1.0, SxiMAbs, Radius(iXiM));
                    for (int l = 0; l < dof; ++l)
                    {
                        LdU[l] += 0.5 * (HdU[l] - lamXiM * mDU2(iXiM)[l]);
                    }
                }
                if (etaM)
                {
                    mModel.JacobianDU(HdU, block, iEtaM, U(iEtaM), mDU2(iEtaM), Seta(i, j - 1, k), -1.0, SetaMAbs, Radius(iEtaM));
                    for (int l = 0; l < dof; ++l)
                    {
                        LdU[l] += 0.5 * (HdU[l] - lamEtaM * mDU2(iEtaM)[l]);
                    }
                }
                if (zetaM)
                {
                    mModel.JacobianDU(HdU, block, iZetaM, U(iZetaM), mDU2(iZetaM), Szeta(i, j, k - 1), -1.0, SzetaMAbs, Radius(iZetaM));
                    for (int l = 0; l < dof; ++l)
                    {
                        LdU[l] += 0.5 * (HdU[l] - lamZetaM * mDU2(iZetaM)[l]);
                    }
                }

                if (!xiM || !etaM || !zetaM)
                {
                    *mD(i, j, k) *= mBDiag;
                }

                if (unsteady)
                {
                    double* uk = U(i, j, k);
                    double* un1 = mU2(i, j, k); // n
                    double* un2 = mU3(i, j, k); // n - 1
                    double* rhs = mRHS(i, j, k);

#if FIRST_ORDER
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] += -vol * (uk[l] - un1[l]) / dTReal;
                    }
                    *mD(i, j, k) += vol / dTReal;
#else
                    for (int l = 0; l < dof; ++l)
                    {
                        rhs[l] += -vol * 0.5 * (3.0 * uk[l] - 4.0 * un1[l] + un2[l]) / dTReal;
                    }
                    *mD(i, j, k) += 1.5 * vol / dTReal;
#endif
                }

                for (int l = 0; l < dof; ++l)
                {
                    mDU2(i, j, k)[l] = (mRHS(i, j, k)[l] - LdU[l]) / (*mD(i, j, k));
                }
            }
        }
    }

    mResidual = Residual(mRHS, block.CellRange(), block.ID());

    delete[] HdU;
    delete[] LdU;
}

template <class Model, class ResEval>
template <class Context>
void
LUSGSIntegrator<Model, ResEval>::Sweep2(const Structured<double>& U)
{
    Structured<double>& dU = DU();

    const Block& block = GetBlock();
    IndexRange cr = block.CellRange();
    const Structured<double>& Sxi   = block.Sxi();
    const Structured<double>& Seta  = block.Seta();
    const Structured<double>& Szeta = block.Szeta();
    const Structured<double>& Radius = block.Radius();

    int dof = mModel.DOF();
    double* HdU = new double[dof];
    double* UdU = new double[dof];

    for (int k = cr.End.K; k >= cr.Start.K; --k)
    {
        for (int j = cr.End.J; j >= cr.Start.J; --j)
        {
            for (int i = cr.End.I; i >= cr.Start.I; --i)
            {
                for (int l = 0; l < dof; ++l)
                {
                    mDU2(i, j, k)[l] *= (*mD(i, j, k));
                }

                bool xiP, etaP, zetaP;
                xiP   = i < cr.End.I;
                etaP  = j < cr.End.J;
                zetaP = k < cr.End.K;

                double SxiPAbs, SetaPAbs, SzetaPAbs;
                SxiPAbs = Vector3(Sxi(i, j, k)).Mag();
                SetaPAbs = Vector3(Seta(i, j, k)).Mag();
                SzetaPAbs = Vector3(Szeta(i, j, k)).Mag();

                IndexIJK i0, iXiP, iEtaP, iZetaP;
                i0 = IndexIJK(i, j, k);
                iXiP   = IndexIJK(i + 1, j, k);
                iEtaP  = IndexIJK(i, j + 1, k);
                iZetaP = IndexIJK(i, j, k + 1);

                double lamXiP, lamEtaP, lamZetaP;
                lamXiP   = mModel.ScalarCoeff(block, U, i0, iXiP,     Sxi(i, j, k), 1.0, SxiPAbs, Radius);
                lamEtaP  = mModel.ScalarCoeff(block, U, i0, iEtaP,   Seta(i, j, k), 1.0, SetaPAbs, Radius);
                lamZetaP = mModel.ScalarCoeff(block, U, i0, iZetaP, Szeta(i, j, k), 1.0, SzetaPAbs, Radius);

                for (int l = 0; l < dof; ++l)
                {
                    UdU[l] = 0.0;
                }
                if (xiP)
                {
                    mModel.JacobianDU(HdU, block, iXiP, U(iXiP), mDU2(iXiP), Sxi(i, j, k), 1.0, SxiPAbs, Radius(iXiP));
                    for (int l = 0; l < dof; ++l)
                    {
                        UdU[l] += 0.5 * (HdU[l] - lamXiP * mDU2(iXiP)[l]);
                    }
                }
                if (etaP)
                {
                    mModel.JacobianDU(HdU, block, iEtaP, U(iEtaP), mDU2(iEtaP), Seta(i, j, k), 1.0, SetaPAbs, Radius(iEtaP));
                    for (int l = 0; l < dof; ++l)
                    {
                        UdU[l] += 0.5 * (HdU[l] - lamEtaP * mDU2(iEtaP)[l]);
                    }
                }
                if (zetaP)
                {
                    mModel.JacobianDU(HdU, block, iZetaP, U(iZetaP), mDU2(iZetaP), Szeta(i, j, k), 1.0, SzetaPAbs, Radius(iZetaP));
                    for (int l = 0; l < dof; ++l)
                    {
                        UdU[l] += 0.5 * (HdU[l] - lamZetaP * mDU2(iZetaP)[l]);
                    }
                }

                for (int l = 0; l < dof; ++l)
                {
                    dU(i, j, k)[l] = (mDU2(i, j, k)[l] - UdU[l]) / (*mD(i, j, k));
                }
            }
        }
    }

    delete[] HdU;
    delete[] UdU;
}

#endif

