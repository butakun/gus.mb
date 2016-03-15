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
// $Id: SolverFunctions.cpp 168 2011-12-13 14:28:38Z kato $

#include "SolverFunctions.h"
#include <cmath>
#include <algorithm>

void UpdateSolution(
    IndexRange cellRange,
    Structured<double> U,
    Structured<double> DT,
    Structured<double> Vol,
    Structured<double> RHS
    )
{
    for (int i = cellRange.Start.I; i < cellRange.End.I; ++i)
    {
        for (int j = cellRange.Start.J; j < cellRange.End.J; ++j)
        {
            for (int k = cellRange.Start.K; k < cellRange.End.K; ++k)
            {
                double* RR = RHS(i, j, k);
                double* UU = U(i, j, k);
                double dt = *DT(i, j, k);
                double vol = *Vol(i, j, k);
                for (int l = 0; l < 5; ++l)
                {
                    UU[l] = UU[l] - RR[l] * dt / vol;
                }
            }
        }
    }
}

void EvaluateFlux(
    double* UL, double* UR,
    double* Sn,
    double gamma,
    double* H
    )
{
    double HL[5], HR[5];
    double sx, sy, sz, sn;
    double rhoL, uL, vL, wL, pL, rhoeL, vnL, rhohtL, cL, lambdaL;
    double rhoR, uR, vR, wR, pR, rhoeR, vnR, rhohtR, cR, lambdaR;
    double lambdaLR;

    sx = Sn[0];
    sy = Sn[1];
    sz = Sn[2];
    sn = std::sqrt(sx * sx + sy * sy  + sz * sz);

    rhoL = UL[0];
    uL = UL[1] / rhoL;
    vL = UL[2] / rhoL;
    wL = UL[3] / rhoL;
    rhoeL = UL[4] - 0.5 * rhoL * (uL * uL + vL * vL + wL * wL);
    pL = (gamma - 1.0) * rhoeL;
    vnL = sx * uL + sy * vL + sz * wL;
    rhohtL = UL[4] + pL;
    cL = std::sqrt(gamma * pL / rhoL);
    lambdaL = std::abs(vnL) + sn * cL;

    rhoR = UR[0];
    uR = UR[1] / rhoR;
    vR = UR[2] / rhoR;
    wR = UR[3] / rhoR;
    rhoeR = UR[4] - 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);
    pR = (gamma - 1.0) * rhoeR;
    vnR = sx * uR + sy * vR + sz * wR;
    rhohtR = UR[4] + pR;
    cR = std::sqrt(gamma * pR / rhoR);
    lambdaR = std::abs(vnR) + sn * cR;

    HL[0] = rhoL * vnL;
    HL[1] = rhoL * uL * vnL + sx * pL;
    HL[2] = rhoL * vL * vnL + sy * pL;
    HL[3] = rhoL * wL * vnL + sz * pL;
    HL[4] = vnL * rhohtL;

    HR[0] = rhoR * vnR;
    HR[1] = rhoR * uR * vnR + sx * pR;
    HR[2] = rhoR * vR * vnR + sy * pR;
    HR[3] = rhoR * wR * vnR + sz * pR;
    HR[4] = vnR * rhohtR;

    H[0] = 0.5 * (HL[0] + HR[0]);
    H[1] = 0.5 * (HL[1] + HR[1]);
    H[2] = 0.5 * (HL[2] + HR[2]);
    H[3] = 0.5 * (HL[3] + HR[3]);
    H[4] = 0.5 * (HL[4] + HR[4]);

    lambdaLR = std::max(lambdaL, lambdaR);
    H[0] -= 0.5 * lambdaLR * (UR[0] - UL[0]);
    H[1] -= 0.5 * lambdaLR * (UR[1] - UL[1]);
    H[2] -= 0.5 * lambdaLR * (UR[2] - UL[2]);
    H[3] -= 0.5 * lambdaLR * (UR[3] - UL[3]);
    H[4] -= 0.5 * lambdaLR * (UR[4] - UL[4]);
}

inline
void ConservativeToPrimitive(double* U, double* Q, double gamma)
{
    // Q = rho, u, v, w, p
    Q[0] = U[0];
    Q[1] = U[1] / U[0];
    Q[2] = U[2] / U[0];
    Q[3] = U[3] / U[0];
    Q[4] = U[4] - 0.5 * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) / U[0];
    Q[4] = Q[4] * (gamma - 1.0);

}

inline
void PrimitiveToConservative(double* Q, double* U, double gamma)
{
    // Q = rho, u, v, w, p
    U[0] = Q[0];
    U[1] = Q[0] * Q[1];
    U[2] = Q[0] * Q[2];
    U[3] = Q[0] * Q[3];
    U[4] = Q[4] / (gamma - 1.0) + 0.5 * Q[0] * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]);
}

inline
double FirstOrder(double r)
{
    return 0.0;
}

inline
double MinMod(double r)
{
    return std::max(0.0, std::min(1.0, r));
}

inline
double VanAlbada(double r)
{
    return r <= 0.0 ? 0.0 : (r * r + r) / (r * r + 1.0);
}

inline
double VanLeer(double r)
{
    return (r + std::abs(r)) / (1.0 + std::abs(r));
}

void Limit(
    double* U1,
    double* U2,
    double* U3,
    double* U4,
    double gamma,
    double* UL,
    double* UR,
    double (*limiter)(double)
    )
{
    double Q1[5], Q2[5], Q3[5], Q4[5];
    ConservativeToPrimitive(U1, Q1, gamma);
    ConservativeToPrimitive(U2, Q2, gamma);
    ConservativeToPrimitive(U3, Q3, gamma);
    ConservativeToPrimitive(U4, Q4, gamma);

    double psiL[5], psiR[5];
    for (int i = 0; i < 5; ++i)
    {
        double rL, rR;
        rL = (Q3[i] - Q2[i]) / (Q2[i] - Q1[i]);
        rR = (Q3[i] - Q2[i]) / (Q4[i] - Q3[i]);
        psiL[i] = limiter(rL);
        psiR[i] = limiter(rR);
    }

    //double psiLMin = std::min(psiL[0], std::min(psiL[1], std::min(psiL[2], std::min(psiL[3], psiL[4]))));
    //double psiRMin = std::min(psiR[0], std::min(psiR[1], std::min(psiR[2], std::min(psiR[3], psiR[4]))));

    double QL[5], QR[5];
    for (int i = 0; i < 5; ++i)
    {
        QL[i] = Q2[i] + 0.5 * psiL[i] * (Q2[i] - Q1[i]);
        QR[i] = Q3[i] - 0.5 * psiR[i] * (Q4[i] - Q3[i]);
    }

    PrimitiveToConservative(QL, UL, gamma);
    PrimitiveToConservative(QR, UR, gamma);
}

void EvaluateResidual(
    double* U,
    double* Sxi, double* Seta,
    int iStart, int iEnd,
    int jStart, int jEnd,
    int idim, int jdim,
    double* R
    )
{
    double H[5];

    for (int i = 0; i < (idim * jdim * 5); ++i)
    {
        R[i] = 0.0;
    }

    for (int i = iStart - 1; i < iEnd; ++i)
    {
        for (int j = jStart; j < jEnd; ++j)
        {
            int ii = i * jdim + j;
            int il = ii, ir  = ii + jdim * 1;
            double *UL, *UR, *Sn;
            double *RL, *RR;
            UL = &U[il * 5];
            UR = &U[ir * 5];
            Sn = &Sxi[ii * 3];
            EvaluateFlux(UL, UR, Sn, 1.4, H);

            RL = &R[il * 5];
            RR = &R[ir * 5];
            for (int l = 0; l < 5; ++l)
            {
                RL[l] += H[l];
                RR[l] -= H[l];
            }
        }
    }

    for (int i = iStart; i < iEnd; ++i)
    {
        for (int j = jStart - 1; j < jEnd; ++j)
        {
            int ii = i * jdim + j;
            int il = ii, ir  = ii + 1;
            double *UL, *UR, *Sn;
            double *RL, *RR;
            UL = &U[il * 5];
            UR = &U[ir * 5];
            Sn = &Seta[ii * 3];
            EvaluateFlux(UL, UR, Sn, 1.4, H);

            RL = &R[il * 5];
            RR = &R[ir * 5];
            for (int l = 0; l < 5; ++l)
            {
                RL[l] += H[l];
                RR[l] -= H[l];
            }
        }
    }
}

void ComputeMetrics(
    double* XYZ,
    double* Sxi, double* Seta, double* Vol,
    int istart, int iend,
    int jstart, int jend,
    int idim, int jdim
    )
{
    // Sxi
    for (int i = istart - 1; i < iend; ++i)
    {
        for (int j = jstart; j < jend; ++j)
        {
            int ii = i * jdim + j;

            int iiS, iiN;
            iiS = i * jdim + j - 1;
            iiN = i * jdim + j;

            double *PS, *PN;
            PS = &XYZ[iiS * 3];
            PN = &XYZ[iiN * 3];

            double *Sn = &Sxi[ii * 3];
            Sn[0] = PN[1] - PS[1];
            Sn[1] = -(PN[0] - PS[0]);
            Sn[2] = 0.0;
        }
    }

    // Seta
    for (int i = istart; i < iend; ++i)
    {
        for (int j = jstart - 1; j < jend; ++j)
        {
            int ii = i * jdim + j;

            int iiW, iiE;
            iiW = (i  - 1) * jdim + j;
            iiE = i * jdim + j;

            double *PW, *PE;
            PW = &XYZ[iiW * 3];
            PE = &XYZ[iiE * 3];

            double *Sn = &Seta[ii * 3];
            Sn[0] = -(PE[1] - PW[1]);
            Sn[1] = PE[0] - PW[0];
            Sn[2] = 0.0;
        }
    }

    // Vol
    for (int i = istart; i < iend; ++i)
    {
        for (int j = jstart; j < jend; ++j)
        {
            int ii = i * jdim + j;
            int iiSW, iiSE, iiNW, iiNE;
            iiSW = (i - 1) * jdim + (j - 1);
            iiSE = (i    ) * jdim + (j - 1);
            iiNW = (i - 1) * jdim + (j    );
            iiNE = (i    ) * jdim + (j    );
            double *PSW, *PSE, *PNW, *PNE;
            PSW = &XYZ[iiSW * 3];
            PSE = &XYZ[iiSE * 3];
            PNW = &XYZ[iiNW * 3];
            PNE = &XYZ[iiNE * 3];
            double d1[3], d2[3];
            d1[0] = PNE[0] - PSW[0];
            d1[1] = PNE[1] - PSW[1];
            d1[2] = PNE[2] - PSW[2];
            d2[0] = PNW[0] - PSE[0];
            d2[1] = PNW[1] - PSE[1];
            d2[2] = PNW[2] - PSE[2];
            Vol[ii] = 0.5 * (d1[0] * d2[1] - d1[1] * d2[0]);
        }
    }
}

void ComputeSpectralRadius(
    double* U,
    double* Sxi, double* Seta, double* Vol,
    double gamma,
    int istart, int iend,
    int jstart, int jend,
    int idim, int jdim,
    double* Lambda
    )
{
    for (int i = istart; i < iend; ++i)
    {
        for (int j = jstart; j < jend; ++j)
        {
            int ii = i * jdim + j;
            int iS = ii - jdim * 1;
            int iW = ii - 1;

            double sxi[3], seta[3];
            sxi[0] = 0.5 * (Sxi[ii * 3 + 0] + Sxi[iS * 3 + 0]);
            sxi[1] = 0.5 * (Sxi[ii * 3 + 1] + Sxi[iS * 3 + 1]);
            sxi[2] = 0.5 * (Sxi[ii * 3 + 2] + Sxi[iS * 3 + 2]);
            seta[0] = 0.5 * (Seta[ii * 3 + 0] + Seta[iW * 3 + 0]);
            seta[1] = 0.5 * (Seta[ii * 3 + 1] + Seta[iW * 3 + 1]);
            seta[2] = 0.5 * (Seta[ii * 3 + 2] + Seta[iW * 3 + 2]);

            double sxiAbs, setaAbs;
            sxiAbs = std::sqrt(sxi[0] * sxi[0] + sxi[1] * sxi[1] + sxi[2] * sxi[2]);
            setaAbs = std::sqrt(seta[0] * seta[0] + seta[1] * seta[1] + seta[2] * seta[2]);

            double* UU = &U[ii * 5];
            double rho, u, v, w, rhoe, p, c;
            rho = UU[0];
            u = UU[1] / rho;
            v = UU[2] / rho;
            w = UU[3] / rho;
            rhoe = UU[4] - 0.5 * rho * (u * u + v * v + w * w);
            p = rhoe * (gamma - 1.0);
            c = std::sqrt(gamma * p / rho);

            double lambdaXi, lambdaEta;
            lambdaXi = std::abs(sxi[0] * u + sxi[1] * v + sxi[2] * w) + sxiAbs * c;
            lambdaEta = std::abs(seta[0] * u + seta[1] * v + seta[2] * w) + setaAbs * c;

            Lambda[ii] = lambdaXi + lambdaEta;
            //Lambda[ii] = std::max(lambdaXi, lambdaEta);
        }
    }
}

void EvaluateResidual(Block& block, Structured<double> R, int order)
{
#if 0
    IndexRange cr = block.CellRange();
    IndexRange dim = block.Dim();
    EvaluateResidual(
        block.U().Data,
        block.Sxi().Data,
        block.Seta().Data,
        cr.Start.I, cr.End.I,
        cr.Start.J, cr.End.J,
        dim.Size().I, dim.Size().J,
        R.Data
        );
#else
    EvaluateResidual(block, block.U(), R, order);
#endif
}

void EvaluateResidual(Block& block, Structured<double> U, Structured<double> R, int order)
{
    IndexRange cr = block.CellRange();
#if 0
    IndexRange dim = block.Dim();
    EvaluateResidual(
        U.Data,
        block.Sxi().Data,
        block.Seta().Data,
        cr.Start.I, cr.End.I,
        cr.Start.J, cr.End.J,
        dim.Size().I, dim.Size().J,
        R.Data
        );
#else
    double (*limiter)(double);
    limiter = MinMod;
    //limiter = VanAlbada;
    //limiter = VanLeer;

    Structured<double> Sxi = block.Sxi();
    Structured<double> Seta = block.Seta();

    double gamma = 1.4;
    double H[5];
    R.SetTo(0.0);

    for (int i = cr.Start.I - 1; i < cr.End.I; ++i)
    {
        for (int j = cr.Start.J; j < cr.End.J; ++j)
        {
            int k = 0;

            double UL[5], UR[5], *Sn;
            double *RL, *RR;
            if (order < 2 || i == cr.Start.I - 1 || i == cr.End.I - 1)
            {
                for (int l = 0; l < 5; ++l)
                {
                    UL[l] = U(i, j, k)[l];
                    UR[l] = U(i + 1, j, k)[l];
                }
            }
            else
            {
                Limit(U(i - 1, j, k), U(i, j, k), U(i + 1, j, k), U(i + 2, j, k), gamma, UL, UR, limiter);
            }
            Sn = Sxi(i, j, k);
            EvaluateFlux(UL, UR, Sn, gamma, H);

            RL = R(i, j, k);
            RR = R(i + 1, j, k);
            for (int l = 0; l < 5; ++l)
            {
                RL[l] += H[l];
                RR[l] -= H[l];
            }
        }
    }

    for (int i = cr.Start.I; i < cr.End.I; ++i)
    {
        for (int j = cr.Start.J - 1; j < cr.End.J; ++j)
        {
            int k = 0;

            double UL[5], UR[5], *Sn;
            double *RL, *RR;
            if (order < 2 || j == cr.Start.J - 1 || j == cr.End.J - 1)
            {
                for (int l = 0; l < 5; ++l)
                {
                    UL[l] = U(i, j, k)[l];
                    UR[l] = U(i, j + 1, k)[l];
                }
            }
            else
            {
                Limit(U(i, j - 1, k), U(i, j, k), U(i, j + 1, k), U(i, j + 2, k), gamma, UL, UR, limiter);
            }
            Sn = Seta(i, j, k);
            EvaluateFlux(UL, UR, Sn, gamma, H);

            RL = R(i, j, k);
            RR = R(i, j + 1, k);
            for (int l = 0; l < 5; ++l)
            {
                RL[l] += H[l];
                RR[l] -= H[l];
            }
        }
    }
#endif
}

