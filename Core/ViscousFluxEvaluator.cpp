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
// $Id: ViscousFluxEvaluator.cpp 168 2011-12-13 14:28:38Z kato $

#include "ViscousFluxEvaluator.h"
#include "Block.h"
#include "Physics.h"
#include "Vector3.h"
#include "AuxiliaryCell.h"
#include "GradientEvaluator.h"

class VelocityTemperatureAdaptor
{
public:
    VelocityTemperatureAdaptor(const Structured<double>& UU_, const Structured<double>& Rad_, double gamma)
    : UU(UU_), Rad(Rad_), Gamma(gamma)
    {
    }

    void Evaluate(double* U, const IndexIJK& i) const
    {
        double* uu = UU(i);
        double* rad = Rad(i);
        double rhoe, p, romegaSq;
        // Q = [u, v, w, T]
        U[0] = uu[1] / uu[0];
        U[1] = uu[2] / uu[0];
        U[2] = uu[3] / uu[0];
        romegaSq = rad[3];
        rhoe = uu[4] - 0.5 * (uu[1] * uu[1] + uu[2] * uu[2] + uu[3] * uu[3]) / uu[0] + 0.5 * uu[0] * romegaSq;
        p = (Gamma - 1.0) * rhoe;
        U[3] = p / uu[0]; // p = rho * T (when using the OVERFLOW nondimensionalization.)
    }

protected:

private:
    const Structured<double>& UU;
    const Structured<double>& Rad;
    double Gamma;
};

ViscousFluxEvaluator::ViscousFluxEvaluator()
{
}

inline
void Velocity(double& u, double& v, double& w, double* U1, double* U2)
{
    u = 0.5 * (U1[1] / U1[0] + U2[1] / U2[0]);
    v = 0.5 * (U1[2] / U1[0] + U2[2] / U2[0]);
    w = 0.5 * (U1[3] / U1[0] + U2[3] / U2[0]);
}

inline
void QuadAverage(double* Q, double* Q1, double* Q2, double* Q3, double* Q4)
{
    Q[0] = 0.25 * (Q1[0] + Q2[0] + Q3[0] + Q4[0]);
    Q[1] = 0.25 * (Q1[1] + Q2[1] + Q3[1] + Q4[1]);
    Q[2] = 0.25 * (Q1[2] + Q2[2] + Q3[2] + Q4[2]);
    Q[3] = 0.25 * (Q1[3] + Q2[3] + Q3[3] + Q4[3]);
}

inline
void ViscousFlux(
    double* Hv, const AuxiliaryCell& cell,
    const Structured<double>& U, const Structured<double>& MuK, const Structured<double>& TurMuK,
    double dUdX[4][3], double gamma, double ReL, double Pr, double PrT
    )
{
    // Cell center indices (the pair sandwiching the cell face of interest)
    IndexIJK iU1M = cell.IndexU1M();
    IndexIJK iU1P = cell.IndexU1P();

    // Properties at cell face
    double mu, mut;
    mu = 0.5 * (MuK(iU1M)[0] + MuK(iU1P)[0]);
    mut = 0.5 * (TurMuK(iU1M)[0] + TurMuK(iU1P)[0]);
    //double k, kt;
    //k  = 0.5 * (MuK(iU1M)[1] + MuK(iU1P)[1]);
    //kt  = 0.5 * (TurMuK(iU1M)[1] + TurMuK(iU1P)[1]);

    double u, v, w;
    Velocity(u, v, w, U(iU1M), U(iU1P));

    double dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dTdx, dTdy, dTdz;
    dudx = dUdX[0][0]; dudy = dUdX[0][1]; dudz = dUdX[0][2];
    dvdx = dUdX[1][0]; dvdy = dUdX[1][1]; dvdz = dUdX[1][2];
    dwdx = dUdX[2][0]; dwdy = dUdX[2][1]; dwdz = dUdX[2][2];
    dTdx = dUdX[3][0]; dTdy = dUdX[3][1]; dTdz = dUdX[3][2];

    double divV;
    divV = dudx + dvdy + dwdz;

    double muu = mu + mut;
    double txx, tyy, tzz, txy, txz, tyz;
    txx = muu / ReL * (2.0 * dudx - 2.0 / 3.0 * divV);
    tyy = muu / ReL * (2.0 * dvdy - 2.0 / 3.0 * divV);
    tzz = muu / ReL * (2.0 * dwdz - 2.0 / 3.0 * divV);
    txy = muu / ReL * (dudy + dvdx);
    txz = muu / ReL * (dudz + dwdx);
    tyz = muu / ReL * (dvdz + dwdy);

    //double qcoeff = mu * gamma / ((gamma - 1.0) * ReL * Pr);
    double qcoeff = gamma / ((gamma - 1.0) * ReL) * (mu / Pr + mut / PrT);
    double qx, qy, qz;
    qx = qcoeff * dTdx;
    qy = qcoeff * dTdy;
    qz = qcoeff * dTdz;

    double phix, phiy, phiz;
    phix = u * (txx + txy + txz) + qx;
    phiy = v * (txy + tyy + tyz) + qy;
    phiz = w * (txz + tyz + tzz) + qz;

    const Vector3& S = cell.Surface();
    double sx = S[0], sy = S[1], sz = S[2];
    Hv[0] = 0.0;
    Hv[1] = sx * txx + sy * txy + sz * txz;
    Hv[2] = sx * txy + sy * tyy + sz * tyz;
    Hv[3] = sx * txz + sy * tyz + sz * tzz;
    Hv[4] = sx * phix + sy * phiy + sz * phiz;
}

void
ViscousFluxEvaluator::EvaluateFlux(
    const Structured<double>& RHS,
    const Block& block,
    const Structured<double>& U
    ) const
{
    IndexRange cr = block.CellRange();

    const Structured<double>& Radius = block.Radius();
    const Structured<double>& MuK = block.MuK();
    const Structured<double>& TurMuK = block.TurMuK();

    double gamma = Physics::GetInstance()->Gamma();
    double ReL = Physics::GetInstance()->ReynoldsNumber();
    double Pr = Physics::GetInstance()->PrandtlNumber();
    double PrT = 0.9; // FIXME

    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I - 1; i <= cr.End.I; ++i)
            {
                bool boundary1 = i == cr.Start.I - 1;
                bool boundary2 = i == cr.End.I;

                IndexIJK iU(i, j, k);

                VelocityTemperatureAdaptor qadaptor(U, Radius, gamma);
                AuxiliaryCell cell(I, iU, block);
                AuxiliaryCellAdaptor<4, VelocityTemperatureAdaptor> adaptor(cell, qadaptor);
                GradientEvaluator<4, AuxiliaryCellAdaptor<4, VelocityTemperatureAdaptor> > gradEval(adaptor);

                // Evaluate dUdX where U = {u, v, w, T} and X = {x, y, z}
                double dUdX[4][3];
                gradEval.Evaluate(dUdX, iU);

                double Hv[5];
                ViscousFlux(Hv, cell, U, MuK, TurMuK, dUdX, gamma, ReL, Pr, PrT);

                if (!boundary1)
                {
                    double* RL = RHS(i, j, k);
                    for (int l = 0; l < 5; ++l)
                    {
                        RL[l] -= Hv[l];
                    }
                }
                if (!boundary2)
                {
                    double* RR = RHS(i + 1, j, k);
                    for (int l = 0; l < 5; ++l)
                    {
                        RR[l] += Hv[l];
                    }
                }
            }
        }
    }

    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J - 1; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                bool boundary1 = j == cr.Start.J - 1;
                bool boundary2 = j == cr.End.J;

                IndexIJK iU(i, j, k);

                VelocityTemperatureAdaptor qadaptor(U, Radius, gamma);
                AuxiliaryCell cell(J, iU, block);
                AuxiliaryCellAdaptor<4, VelocityTemperatureAdaptor> adaptor(cell, qadaptor);
                GradientEvaluator<4, AuxiliaryCellAdaptor<4, VelocityTemperatureAdaptor> > gradEval(adaptor);

                // Evaluate dUdX where U = {u, v, w, T} and X = {x, y, z}
                double dUdX[4][3];
                gradEval.Evaluate(dUdX, iU);

                double Hv[5];
                ViscousFlux(Hv, cell, U, MuK, TurMuK, dUdX, gamma, ReL, Pr, PrT);

                if (!boundary1)
                {
                    double* RL = RHS(i, j, k);
                    for (int l = 0; l < 5; ++l)
                    {
                        RL[l] -= Hv[l];
                    }
                }
                if (!boundary2)
                {
                    double* RR = RHS(i, j + 1, k);
                    for (int l = 0; l < 5; ++l)
                    {
                        RR[l] += Hv[l];
                    }
                }
            }
        }
    }

    for (int k = cr.Start.K - 1; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                bool boundary1 = k == cr.Start.K - 1;
                bool boundary2 = k == cr.End.K;

                IndexIJK iU(i, j, k);

                VelocityTemperatureAdaptor qadaptor(U, Radius, gamma);
                AuxiliaryCell cell(K, iU, block);
                AuxiliaryCellAdaptor<4, VelocityTemperatureAdaptor> adaptor(cell, qadaptor);
                GradientEvaluator<4, AuxiliaryCellAdaptor<4, VelocityTemperatureAdaptor> > gradEval(adaptor);

                // Evaluate dUdX where U = {u, v, w, T} and X = {x, y, z}
                double dUdX[4][3];
                gradEval.Evaluate(dUdX, iU);

                double Hv[5];
                ViscousFlux(Hv, cell, U, MuK, TurMuK, dUdX, gamma, ReL, Pr, PrT);

                if (!boundary1)
                {
                    double* RL = RHS(i, j, k);
                    for (int l = 0; l < 5; ++l)
                    {
                        RL[l] -= Hv[l];
                    }
                }
                if (!boundary2)
                {
                    double* RR = RHS(i, j, k + 1);
                    for (int l = 0; l < 5; ++l)
                    {
                        RR[l] += Hv[l];
                    }
                }
            }
        }
    }
}

