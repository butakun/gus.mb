// $Id: BCViscousWall.cpp 277 2013-06-04 01:58:51Z kato $

#include "BCViscousWall.h"
#include "Physics.h"
#include "Vector3.h"

BCViscousWall::BCViscousWall(const IndexRange& meshRange, Direction dir, double TWall)
:   BCPlanarLocal(meshRange, dir),
    mTWall(TWall / Physics::GetInstance()->TRef())
{
    if (TWall < 0.0)
        mTWall = -1.0;
}

void
BCViscousWall::SetMask(Structured<int>& mask) const
{
    SetMaskWithAValue(mask, 1);
}

void
BCViscousWall::LocalFunc(
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

    double gamma, pi, Ti, romegaSqi, rhokeri;
    gamma = Physics::GetInstance()->Gamma();
    romegaSqi = Ri[3];
    rhokeri = 0.5 * Ui[0] * romegaSqi;
    pi = (gamma - 1.0) * (Ui[4] - 0.5 * (Ui[1] * Ui[1] + Ui[2] * Ui[2] + Ui[3] * Ui[3]) / Ui[0] + rhokeri);
    Ti = pi / Ui[0];

    double rhoGhost, uGhost, vGhost, wGhost, TGhost, pGhost, rhoetGhost;
    double romegaSqg, rhokerg;
    rhoGhost = Ui[0];
    uGhost = -Ui[1] / Ui[0];
    vGhost = -Ui[2] / Ui[0];
    wGhost = -Ui[3] / Ui[0];
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
BCViscousWall::LocalFuncTurb(
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

