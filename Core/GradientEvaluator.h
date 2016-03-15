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
// $Id: GradientEvaluator.h 55 2010-07-27 14:21:56Z kato $
#ifndef INCLUDED_GRADIENT_EVALUATOR_H__
#define INCLUDED_GRADIENT_EVALUATOR_H__

template <int DOF, class QADAPTOR>
class PrimaryCellAdaptor
{
public:
    PrimaryCellAdaptor(
        const QADAPTOR& qAdaptor,
        const Block& block
        )
    : DI(1, 0, 0), DJ(0, 1, 0), DK(0, 0, 1), QAdaptor(qAdaptor), mBlock(block)
    {}

    void SetIndex(const IndexIJK& ijk) { I = ijk; }
    Vector3 Surface1M() const { return Vector3(mBlock.Sxi()(I - DI)); }
    Vector3 Surface1P() const { return Vector3(mBlock.Sxi()(I)); }
    Vector3 Surface2M() const { return Vector3(mBlock.Seta()(I - DJ)); }
    Vector3 Surface2P() const { return Vector3(mBlock.Seta()(I)); }
    Vector3 Surface3M() const { return Vector3(mBlock.Szeta()(I - DK)); }
    Vector3 Surface3P() const { return Vector3(mBlock.Szeta()(I)); }
    double Volume() const { return *mBlock.Vol()(I); }
    void U1M(double* U) const { Average(U, I, I - DI); }
    void U1P(double* U) const { Average(U, I, I + DI); }
    void U2M(double* U) const { Average(U, I, I - DJ); }
    void U2P(double* U) const { Average(U, I, I + DJ); }
    void U3M(double* U) const { Average(U, I, I - DK); }
    void U3P(double* U) const { Average(U, I, I + DK); }

protected:
    void Average(double* U, const IndexIJK& i1, const IndexIJK& i2) const
    {
        double U1[DOF], U2[DOF];
        QAdaptor.Evaluate(U1, i1);
        QAdaptor.Evaluate(U2, i2);
        for (int l = 0; l < DOF; ++l)
        {
            U[l] = 0.5 * (U1[l] + U2[l]);
        }
    }

private:
    const IndexIJK DI, DJ, DK;

    QADAPTOR QAdaptor;
    IndexIJK I;
    const Block& mBlock;
};

template <int DOF, class QADAPTOR>
class AuxiliaryCellAdaptor
{
public:
    AuxiliaryCellAdaptor(const AuxiliaryCell& cell, const QADAPTOR& qAdaptor)
    : Cell(cell), QAdaptor(qAdaptor)
    {}

    void SetIndex(const IndexIJK& ijk) { I = ijk; Cell.SetIndex(ijk); }
    const Vector3& Surface1M() const { return Cell.Surface1M(); }
    const Vector3& Surface1P() const { return Cell.Surface1P(); }
    const Vector3& Surface2M() const { return Cell.Surface2M(); }
    const Vector3& Surface2P() const { return Cell.Surface2P(); }
    const Vector3& Surface3M() const { return Cell.Surface3M(); }
    const Vector3& Surface3P() const { return Cell.Surface3P(); }
    double Volume() const { return Cell.Volume(); }
    void U1M(double* U) const { QAdaptor.Evaluate(U, Cell.IndexU1M()); }
    void U1P(double* U) const { QAdaptor.Evaluate(U, Cell.IndexU1P()); }
    void U2M(double* U) const { IndexIJK i1, i2, i3, i4; Cell.IndicesU2M(i1, i2, i3, i4); QuadAverage(U, i1, i2, i3, i4); }
    void U2P(double* U) const { IndexIJK i1, i2, i3, i4; Cell.IndicesU2P(i1, i2, i3, i4); QuadAverage(U, i1, i2, i3, i4); }
    void U3M(double* U) const { IndexIJK i1, i2, i3, i4; Cell.IndicesU3M(i1, i2, i3, i4); QuadAverage(U, i1, i2, i3, i4); }
    void U3P(double* U) const { IndexIJK i1, i2, i3, i4; Cell.IndicesU3P(i1, i2, i3, i4); QuadAverage(U, i1, i2, i3, i4); }

protected:
    void QuadAverage(double* U, const IndexIJK& i1, const IndexIJK& i2, const IndexIJK& i3, const IndexIJK& i4) const
    {
        double U1[DOF], U2[DOF], U3[DOF], U4[DOF];
        QAdaptor.Evaluate(U1, i1);
        QAdaptor.Evaluate(U2, i2);
        QAdaptor.Evaluate(U3, i3);
        QAdaptor.Evaluate(U4, i4);
        for (int l = 0; l < DOF; ++l)
        {
            U[l] = 0.25 * (U1[l] + U2[l] + U3[l] + U4[l]);
        }
    }

private:
    IndexIJK I;
    AuxiliaryCell Cell;
    QADAPTOR QAdaptor;
};

template <int DOF, class ADAPTOR>
class GradientEvaluator
{
public:
    GradientEvaluator(const ADAPTOR& adaptor);

    void Evaluate(double dUdX[DOF][3], const IndexIJK& ijk);

protected:

private:
    ADAPTOR Adaptor;
};

#include "Vector3.h"

template <int DOF, class ADAPTOR>
inline
GradientEvaluator<DOF, ADAPTOR>::GradientEvaluator(const ADAPTOR& adaptor)
:   Adaptor(adaptor)
{
}

template <int DOF, class ADAPTOR>
inline
void
GradientEvaluator<DOF, ADAPTOR>::Evaluate(double dUdX[][3], const IndexIJK& ijk)
{
    Adaptor.SetIndex(ijk);

    Vector3 S1M = Adaptor.Surface1M();
    Vector3 S1P = Adaptor.Surface1P();
    Vector3 S2M = Adaptor.Surface2M();
    Vector3 S2P = Adaptor.Surface2P();
    Vector3 S3M = Adaptor.Surface3M();
    Vector3 S3P = Adaptor.Surface3P();
    double Vol = Adaptor.Volume();

    double U1M[DOF], U1P[DOF], U2M[DOF], U2P[DOF], U3M[DOF], U3P[DOF];
    Adaptor.U1M(U1M);
    Adaptor.U1P(U1P);
    Adaptor.U2M(U2M);
    Adaptor.U2P(U2P);
    Adaptor.U3M(U3M);
    Adaptor.U3P(U3P);

    for (int l = 0; l < DOF; ++l)
    {
        dUdX[l][0]  = S1P[0] * U1P[l] / Vol;
        dUdX[l][1]  = S1P[1] * U1P[l] / Vol;
        dUdX[l][2]  = S1P[2] * U1P[l] / Vol;
        dUdX[l][0] -= S1M[0] * U1M[l] / Vol;
        dUdX[l][1] -= S1M[1] * U1M[l] / Vol;
        dUdX[l][2] -= S1M[2] * U1M[l] / Vol;
        dUdX[l][0] += S2P[0] * U2P[l] / Vol;
        dUdX[l][1] += S2P[1] * U2P[l] / Vol;
        dUdX[l][2] += S2P[2] * U2P[l] / Vol;
        dUdX[l][0] -= S2M[0] * U2M[l] / Vol;
        dUdX[l][1] -= S2M[1] * U2M[l] / Vol;
        dUdX[l][2] -= S2M[2] * U2M[l] / Vol;
        dUdX[l][0] += S3P[0] * U3P[l] / Vol;
        dUdX[l][1] += S3P[1] * U3P[l] / Vol;
        dUdX[l][2] += S3P[2] * U3P[l] / Vol;
        dUdX[l][0] -= S3M[0] * U3M[l] / Vol;
        dUdX[l][1] -= S3M[1] * U3M[l] / Vol;
        dUdX[l][2] -= S3M[2] * U3M[l] / Vol;
    }
}

#endif // INCLUDED_GRADIENT_EVALUATOR_H__
