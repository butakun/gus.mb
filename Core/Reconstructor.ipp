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
// $Id: Reconstructor.ipp 266 2013-01-31 07:23:08Z kato $

#include "Physics.h"

inline
PrimitiveVariableCodec::PrimitiveVariableCodec()
:   mGamma(Physics::GetInstance()->Gamma())
{
}

inline
void
PrimitiveVariableCodec::Decode(double* Q, double* U, double* Rad) const
{
    double rho, u, v, w, rhoet, rhoke, rhoker, rhoe, p, romegaSq;
    rho = U[0];
    u = U[1] / U[0];
    v = U[2] / U[0];
    w = U[3] / U[0];
    rhoet = U[4];
    romegaSq = Rad[3];

    rhoke = 0.5 * rho * (u * u + v * v + w * w);
    rhoker = 0.5 * rho * romegaSq;
    rhoe = rhoet - rhoke + rhoker;
    p = (mGamma - 1.0) * rhoe;

    Q[0] = rho;
    Q[1] = u;
    Q[2] = v;
    Q[3] = w;
    Q[4] = p;
}

inline
void
PrimitiveVariableCodec::Encode(double* U, double* Q, double* Rad) const
{
    double rho, u, v, w, p, rhoe, rhoke, rhoker, rhoet, romegaSq;
    rho = Q[0];
    u = Q[1];
    v = Q[2];
    w = Q[3];
    p = Q[4];
    romegaSq = Rad[3];
    rhoe = p / (mGamma - 1.0); 
    rhoke = 0.5 * rho * (u * u + v * v + w * w);
    rhoker = 0.5 * rho * romegaSq;
    rhoet = rhoe + rhoke - rhoker;

    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = rho * w;
    U[4] = rhoet;
}

template <int N, class LIMITER, class CODEC>
inline
Reconstructor<N, LIMITER, CODEC>::Reconstructor()
{
}

template <int N, class LIMITER, class CODEC>
inline
void
Reconstructor<N, LIMITER, CODEC>::Reconstruct(
        double* UL, double* UR,
        const Structured<double>& U,
        const Structured<double>& Rad,
        const IndexIJK& I, const IndexIJK& dI
        ) const
{
    double Q1[N], Q2[N], Q3[N], Q4[N], QL[N], QR[N];
    mCodec.Decode(Q1, U(I - dI), Rad(I - dI));
    mCodec.Decode(Q2, U(I), Rad(I));
    mCodec.Decode(Q3, U(I + dI), Rad(I + dI));
    mCodec.Decode(Q4, U(I + dI + dI), Rad(I + dI + dI));

    for (int i = 0; i < N; ++i)
    {
        double d12, d23, d34, rL, rR;
	d12 = Q2[i] - Q1[i];
	d23 = Q3[i] - Q2[i];
	d34 = Q4[i] - Q3[i];
	rL = d23 / d12;
	rR = d23 / d34;

        double phiL, phiR;
        phiL = mLimiter(rL);
        phiR = mLimiter(rR);

        QL[i] = Q2[i] + 0.5 * phiL * (Q2[i] - Q1[i]);
        QR[i] = Q3[i] - 0.5 * phiR * (Q4[i] - Q3[i]);
    }

    // FIXME: rad
    mCodec.Encode(UL, QL, Rad(I));
    mCodec.Encode(UR, QR, Rad(I + dI));
}

template <int N>
inline
FirstOrderReconstructor<N>::FirstOrderReconstructor()
{
}

template <int N>
inline
void
FirstOrderReconstructor<N>::Reconstruct(
        double* UL, double* UR,
        const Structured<double>& U,
        const Structured<double>& Rad,
        const IndexIJK& I, const IndexIJK& dI
        ) const
{
    double* ul = U(I);
    double* ur = U(I + dI);
    for (int i = 0; i < N; ++i)
    {
        UL[i] = ul[i];
        UR[i] = ur[i];
    }
}

