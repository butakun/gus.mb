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
#include "BCRotatingViscousWall.h"
#include "Physics.h"
#include "Vector3.h"

BCRotatingViscousWall::BCRotatingViscousWall(const IndexRange& meshRange, Direction dir, double TWall, const Vector3& angularVelocity)
:   BCPlanarLocal(meshRange, dir),
    mTWall(TWall / Physics::GetInstance()->TRef()),
    mRBM(new RotationalMotion(Vector3(0.0, 0.0, 0.0), angularVelocity))
{
    if (TWall < 0.0)
        mTWall = -1.0;
}

BCRotatingViscousWall::~BCRotatingViscousWall()
{
    delete mRBM;
}

void
BCRotatingViscousWall::SetMask(Structured<int>& mask) const
{
    SetMaskWithAValue(mask, 1);
}

void
BCRotatingViscousWall::LocalFunc(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& U, const Block& block
    )
{
    IndexIJK iI = iInterior;
    IndexIJK iG = iGhost;
    double* Ui = U(iI);
    double* Ri = block.Radius()(iInterior);
    double* Rg = block.Radius()(iGhost);

    Vector3 posWall = 0.5 * (Vector3(Rg[0], Rg[1], Rg[2]) + Vector3(Ri[0], Ri[1], Ri[2]));
    Vector3 velWall = mRBM->VelocityAt(posWall);

    double gamma, pi, Ti, romegaSqi, rhokeri;
    gamma = Physics::GetInstance()->Gamma();
    romegaSqi = Ri[3];
    rhokeri = 0.5 * Ui[0] * romegaSqi;
    pi = (gamma - 1.0) * (Ui[4] - 0.5 * (Ui[1] * Ui[1] + Ui[2] * Ui[2] + Ui[3] * Ui[3]) / Ui[0] + rhokeri);
    Ti = pi / Ui[0];

    double rhoGhost, uGhost, vGhost, wGhost, TGhost, pGhost, rhoetGhost;
    double romegaSqg, rhokerg;
    rhoGhost = Ui[0];
    uGhost = 2.0 * velWall.X() - Ui[1] / Ui[0];
    vGhost = 2.0 * velWall.Y() - Ui[2] / Ui[0];
    wGhost = 2.0 * velWall.Z() - Ui[3] / Ui[0];
    if (mTWall < 0.0)
        TGhost = Ti;
    else
        TGhost = std::max(0.5 * mTWall, 2.0 * mTWall - Ti);
    pGhost = rhoGhost * TGhost;
    romegaSqg = Rg[3];
    rhokerg = 0.5 * rhoGhost * romegaSqg;
    rhoetGhost = pGhost / (gamma - 1.0)
        + 0.5 * rhoGhost * (uGhost * uGhost + vGhost * vGhost + wGhost * wGhost) - rhokerg;

    double UG[5];
    UG[0] = rhoGhost;
    UG[1] = rhoGhost * uGhost;
    UG[2] = rhoGhost * vGhost;
    UG[3] = rhoGhost * wGhost;
    UG[4] = rhoetGhost;

    double dUiG[5]; // UG - Ui, so that UG = Ui + dUiG, UG2 = UG + dUiG, where UG2 is the second rind layer cell.
    dUiG[0] = UG[0] - Ui[0];
    dUiG[1] = UG[1] - Ui[1];
    dUiG[2] = UG[2] - Ui[2];
    dUiG[3] = UG[3] - Ui[3];
    dUiG[4] = UG[4] - Ui[4];

    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* UGhost = U(iG);

        for (int l = 0; l < 5; l++)
        {
            UGhost[l] = Ui[l] + double(i + 1) * dUiG[l];
        }
        iG -= deltaInterior;
    }
}

void
BCRotatingViscousWall::LocalFuncTurb(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& UT, const Structured<double>& U, const Block& block
    )
{
    //mModel.BCWall(iFace, dGhost, dInterior, Sn, UT, block);

    IndexIJK iG = iGhost, iI = iInterior;
    double volInterior = block.Vol()(iI)[0];
    double dy = volInterior / Sn.Mag();

    double mu = 0.5 * (block.MuK()(iG)[0] + block.MuK()(iI)[0]);
    //double rho = 0.5 * (block.U()(iG)[0] + block.U()(iI)[0]);
    double* UTi = UT(iI);

    const double BETA1 = 0.0750;
    double Re = Physics::GetInstance()->ReynoldsNumber();
    double rhoomega = 10.0 * 6.0 * mu / (BETA1 * dy * dy) / (Re * Re); // From Menter, AIAAJ Vol. 32, No. 8, 1994.
    double rhotke = 0.0;

    double UTG[2];

    UTG[0] = 2.0 * rhotke - UTi[0];
    UTG[1] = 2.0 * rhoomega - UTi[1];

    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* UTGhost = UT(iG);

        for (int l = 0; l < 2; l++)
        {
            UTGhost[l] = UTG[l];
        }
        iG -= deltaInterior;
    }
}

