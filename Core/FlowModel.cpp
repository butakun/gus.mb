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
// $Id: FlowModel.cpp 306 2013-10-02 07:03:25Z kato $

#include "Roster.h"
#include "FlowModel.h"
#include "Physics.h"
#include "SolverUtils.h"
#include <cmath>

void
FlowModel::Initialize()
{
    Roster::GetInstance()->CreateNewTag("U");
    Roster::GetInstance()->CreateNewTag("UT");
}

FlowModel::FlowModel()
:   Gamma(Physics::GetInstance()->Gamma()),
    mHartenEps(0.0) // FIXME
{
}

void
FlowModel::ApplyBCs(Block& block) const
{
    block.ApplyBCs();
}

double
FlowModel::ComputeTimeStep(
    const Block& block,
    const Structured<double>& U,
    Structured<double>& DT,
    double cfl,
    bool localTimeStepping
    )
{
    return SolverUtils::ComputeTimeStep(block, U, DT, cfl, localTimeStepping);
}

double
FlowModel::ScalarCoeff(
    const Block& block, const Structured<double>& U,
    const IndexIJK& Ii, const IndexIJK& Ij,
    double* Sij, double SijSign, double SijAbs, const Structured<double>& radius
    ) const
{
    double lambdaI, lambdaV = 0.0;
    lambdaI = ScalarCoeff(U(Ii), U(Ij), Sij, SijSign, SijAbs, radius(Ii), radius(Ij));

    double rhoij, muij, Re;
    rhoij = 0.5 * (U(Ii)[0] + U(Ij)[0]);
    muij = 0.5 * (block.MuK()(Ii)[0] + block.TurMuK()(Ii)[0] + block.MuK()(Ij)[0] + block.TurMuK()(Ij)[0]);
    Re = Physics::GetInstance()->ReynoldsNumber();

    double vol, dx;
    vol = 0.5 * (*block.Vol()(Ii) + *block.Vol()(Ij));
    dx = vol / SijAbs;

    lambdaV = 1.0 * SijAbs * muij / rhoij / Re / dx; // FIXME

    return lambdaI + lambdaV;
}

// Hij = 0.5 * (Hi + Hj) - 0.5 * lambda * (Uj - Ui)
// In structured mesh, Hi = HL, Hj = HR
double
FlowModel::ScalarCoeff(double* Ui, double* Uj, double* Sij, double SijSign, double SijAbs, double* Ri, double* Rj) const
{
    //double SijAbs;
    //SijAbs = std::sqrt(Sij[0] * Sij[0] + Sij[1] * Sij[1] + Sij[2] * Sij[2]);

    double rhoi, rhoj, rhokeri, rhokerj, rhokei, rhokej, rhoei, rhoej;
    double pi, pj, ci, cj, cij, Cij, Vi, Vj, Vij;
    rhoi = Ui[0];
    rhoj = Uj[0];
    rhokeri = 0.5 * rhoi * Ri[3];
    rhokerj = 0.5 * rhoj * Rj[3];
    rhokei  = 0.5 * (Ui[1] * Ui[1] + Ui[2] * Ui[2] + Ui[3] * Ui[3]) / rhoi;
    rhokej  = 0.5 * (Uj[1] * Uj[1] + Uj[2] * Uj[2] + Uj[3] * Uj[3]) / rhoj;
    rhoei = Ui[4] - rhokei + rhokeri;
    rhoej = Uj[4] - rhokej + rhokerj;
    pi = (Gamma - 1.0) * rhoei;
    pj = (Gamma - 1.0) * rhoej;

    ci = std::sqrt(Gamma * pi / rhoi);
    cj = std::sqrt(Gamma * pj / rhoj);
    cij = 0.5 * (ci + cj);
    Cij = SijAbs * cij;

    // SijSign has no significance here. But this parameter is kept in the interface nonetheless,
    // some other model might need to depend on the direction of the flux surface normal.
    Vi = SijSign * (Sij[0] * Ui[1] + Sij[1] * Ui[2] + Sij[2] * Ui[3]) / rhoi;
    Vj = SijSign * (Sij[0] * Uj[1] + Sij[1] * Uj[2] + Sij[2] * Uj[3]) / rhoj;
    Vij = 0.5 * (Vi + Vj);

    double lambdaAbs;
    lambdaAbs = std::abs(Vij) + Cij;
    lambdaAbs = std::max(lambdaAbs, SijAbs * mHartenEps);

    return lambdaAbs;
}

void Flux(double* H, double* U, double* Sn, double SnSign, double SnAbs, double* Radius, double gamma)
{
    double rho, u, v, w, Vn, p, romegaSq, ker;
    rho = U[0];
    u = U[1] / rho;
    v = U[2] / rho;
    w = U[3] / rho;
    Vn = SnSign * (Sn[0] * u + Sn[1] * v + Sn[2] * w);
    romegaSq = Radius[3];
    ker = 0.5 * romegaSq;
    p = (gamma - 1.0) * (U[4] - 0.5 * rho * (u * u + v * v + w * w) + 0.5 * rho * ker);

    H[0] = rho * Vn;
    H[1] = rho * Vn * u + SnAbs * p;
    H[2] = rho * Vn * v + SnAbs * p;
    H[3] = rho * Vn * w + SnAbs * p;
    H[4] = (U[4] + p) * Vn;
}

void
FlowModel::JacobianDU(
    double* JacDU, const Block& block, const IndexIJK& i,
    double* U, double* dU, double* Sij, double SijSign, double SijAbs, double* Radius
    ) const
{
    double eps = 1.0e-8;
    double H[5], H2[5], U2[5];

    for (int l = 0; l < 5; ++l)
        U2[l] = U[l] + eps * dU[l];

    Flux(H, U, Sij, SijSign, SijAbs, Radius, Gamma);
    Flux(H2, U2, Sij, SijSign, SijAbs, Radius, Gamma);

    for (int l = 0; l < 5; ++l)
        JacDU[l] = (H2[l] - H[l]) / eps;
}

void
FlowModel::UplusDU(Structured<double>& U, const Structured<double>& dU, const IndexRange& cellRange, const Block& block) const
{
    for (int k = cellRange.Start.K; k <= cellRange.End.K; ++k)
    {
        for (int j = cellRange.Start.J; j <= cellRange.End.J; ++j)
        {
            for (int i = cellRange.Start.I; i <= cellRange.End.I; ++i)
            {
                double* u = U(i, j, k);
                double* du = dU(i, j, k);

                double rho, rhou, rhov, rhow, rhoet;
                rho = u[0] + du[0];
                rhou = u[1] + du[1];
                rhov = u[2] + du[2];
                rhow = u[3] + du[3];
                rhoet = u[4] + du[4];
#if 0
                double rhoke, rhoe;
                rhoke = 0.5 * (rhou * rhou + rhov * rhov + rhow * rhow) / rho;
                rhoe = rhoet - rhoke; // FIXME: rothalpy
                if (rhoe < 0.0)
                {
                    LOG << "UplusDU:" << IndexIJK(i, j, k) << " rhoe = " << rhoe << " rhoke = " << rhoke << " limiting kinetic energy" << std::endl;
                    rhoet += 1.1 * (-rhoe);
                }
#endif
                u[0] = rho;
                u[1] = rhou;
                u[2] = rhov;
                u[3] = rhow;
                u[4] = rhoet;
            }
        }
    }
}

void
FlowModel::FromGlobalToLocal(double* ULocal, double* UGlobal, const Block& block, const IndexIJK& ijk) const
{
    const RotationalMotion* rm = dynamic_cast<const RotationalMotion*>(block.GetRigidBodyMotion());
    if (rm == NULL)
    {
        for (int l = 0; l < 5; ++l)
            ULocal[l] = UGlobal[l];
        return;
    }

    const Structured<double>& Rad = block.Radius();
    double* rad = Rad(ijk);
    Vector3 pC(rad[0], rad[1], rad[2]);
    double romegaSq = rad[3];

    Vector3 ve = rm->EntrainmentVelocityAt(pC);

    double rho = UGlobal[0];
    Vector3 vg(UGlobal[1] / rho, UGlobal[2] / rho, UGlobal[3] / rho);
    Vector3 vl = vg - ve;
    double rhoE = UGlobal[4] - 0.5 * rho * vg.MagSq();
    double rhoEtLocal = rhoE + 0.5 * rho * (vl.MagSq() - romegaSq);

    ULocal[0] = rho;
    ULocal[1] = rho * vl.X();
    ULocal[2] = rho * vl.Y();
    ULocal[3] = rho * vl.Z();
    ULocal[4] = rhoEtLocal;
}

void
FlowModel::FromLocalToGlobal(double* UGlobal, double* ULocal, const Block& block, const IndexIJK& ijk) const
{
    const RotationalMotion* rm = dynamic_cast<const RotationalMotion*>(block.GetRigidBodyMotion());
    if (rm == NULL)
    {
        for (int l = 0; l < 5; ++l)
            UGlobal[l] = ULocal[l];
        return;
    }

    const Structured<double>& Rad = block.Radius();
    double* rad = Rad(ijk);
    Vector3 pC(rad[0], rad[1], rad[2]);
    double romegaSq = rad[3];

    Vector3 ve = rm->EntrainmentVelocityAt(pC);

    double rho = ULocal[0];
    Vector3 vl(ULocal[1] / rho, ULocal[2] / rho, ULocal[3] / rho);
    Vector3 vg = vl + ve;

    double rhoe = ULocal[4] - 0.5 * rho * vl.MagSq() + 0.5 * romegaSq;
    double rhoetGlobal = rhoe + 0.5 * rho * vg.MagSq();

    UGlobal[0] = rho;
    UGlobal[1] = rho * vg.X();
    UGlobal[2] = rho * vg.Y();
    UGlobal[3] = rho * vg.Z();
    UGlobal[4] = rhoetGlobal;
}

