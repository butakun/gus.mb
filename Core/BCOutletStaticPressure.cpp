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
// $Id: BCOutletStaticPressure.cpp 277 2013-06-04 01:58:51Z kato $

#include "Communicator.h"
#include "BCOutletStaticPressure.h"
#include "Physics.h"

BCOutletStaticPressure::BCOutletStaticPressure(const IndexRange& meshRange, Direction direction, double pressure)
:   BCPlanarLocal(meshRange, direction),
    mStaticPressure(pressure / Physics::GetInstance()->PRef())
{
    Communicator::GetInstance()->Console() << "BCOutletStaticPressure: p(nondim) = " << mStaticPressure << ", PRef = " << Physics::GetInstance()->PRef() << std::endl;
    std::cout << "BCOutletStaticPressure: p(nondim) = " << mStaticPressure << ", PRef = " << Physics::GetInstance()->PRef() << std::endl;
}

void
BCOutletStaticPressure::LocalFunc(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& U, const Block& block
    )
{
    double* Ui = U(iInterior);

    double gamma;
    gamma = Physics::GetInstance()->Gamma();

    double* Rg = block.Radius()(iGhost);
    double* Ri = block.Radius()(iInterior);

    // Ghost cell quantities
    double rhoGhost, uGhost, vGhost, wGhost, pGhost, rhoke, rhoetGhost;
    double romegaSq, rhoker;
    rhoGhost = Ui[0];
    uGhost = Ui[1] / Ui[0];
    vGhost = Ui[2] / Ui[0];
    wGhost = Ui[3] / Ui[0];
    romegaSq = Rg[3];

    // Interior cell quantities
    double rhoi, ui, vi, wi, romegaSqi, pi, vni;
    rhoi = Ui[0];
    ui = Ui[1] / Ui[0];
    vi = Ui[2] / Ui[0];
    wi = Ui[3] / Ui[0];
    romegaSqi = Ri[3];
    pi = (gamma - 1.0) * (Ui[4] - 0.5 * rhoi * (ui * ui + vi * vi + wi * wi) + 0.5 * rhoi * romegaSqi);
    vni = Sn.X() * ui + Sn.Y() * vi + Sn.Z() * wi; // velocity component normal to boundary face, positive means inward (ghost to interior).
    if (vni > 0.0)
    {
        // this indicates flow reversal, so we apply lower static pressure on the ghost cell
        pGhost = 0.9 * pi;
    }
    else
    {
        pGhost = std::min(1.1 * pi, mStaticPressure);
    }

    rhoker = 0.5 * rhoGhost * romegaSq;
    rhoke = 0.5 * rhoGhost * (uGhost * uGhost + vGhost * vGhost + wGhost * wGhost);
    rhoetGhost = pGhost / (gamma - 1.0) + rhoke - rhoker;

    double UG[5];
    UG[0] = rhoGhost;
    UG[1] = rhoGhost * uGhost;
    UG[2] = rhoGhost * vGhost;
    UG[3] = rhoGhost * wGhost;
    UG[4] = rhoetGhost;

    IndexIJK iG = iGhost;
    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* UGhost = U(iG);

        for (int l = 0; l < 5; l++)
        {
            UGhost[l] = UG[l];
        }
        iG -= deltaInterior;
    }
}

void
BCOutletStaticPressure::LocalFuncTurb(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& UT, const Structured<double>& U, const Block& block
    )
{
    IndexIJK iG = iGhost;

    double* UTI = UT(iInterior);
    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* UTG = UT(iG);
        for (int l = 0; l < UT.DOF(); ++l)
        {
            UTG[l] = UTI[l];
        }
        iG -= deltaInterior;
    }
}
