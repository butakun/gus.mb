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

#include "Communicator.h"
#include "BCInletTotalMassFlow.h"
#include "Physics.h"
#include "TurbulenceSpec.h"
#include "RigidBodyMotion.h"
#include <cassert>

BCInletTotalMassFlow::BCInletTotalMassFlow(
    const IndexRange& meshRange, Direction direction,
    double pt, double tt, double mdotFamily,
    int ndoft, const TurbulenceSpec* turbSpec
    )
:   BCPlanarLocal(meshRange, direction),
    mTotalPressure(pt / Physics::GetInstance()->PRef()),
    mTotalTemperature(tt / Physics::GetInstance()->TRef()),
    mMassFlowRateFamilyTarget(-1.0),
    mMassFlowRateFamilyCurrent(-1.0),
    mMassFlowRateLocal(-1.0),
    mTurbFix(new double[ndoft])
{
    assert(ndoft == 2);

    // FIXME: how do we make this model-neutral?
    mTurbFix[0] = turbSpec->TKE_Nondimensional();
    mTurbFix[1] = turbSpec->Omega_Nondimensional();

    Physics* phys = Physice::GetInstance();
    mMassFlowRateFamilyTarget = mdotFamily
        / (phys->RhoRef() * phys->VRef() * phys->LRef() * phys->LRef());

    Communicator::GetInstance()->Console()
        << "BCInletTotalMassFlow: pt(nondim) = " << mTotalPressure
        << ", Tt(nondim) = " << mTotalTemperature
        << ", MassFlowRateFamily = " << mMassFlowRateFamily << std::endl;
}

BCInletTotalMassFlow::~BCInletTotalMassFlow()
{
    delete[] mTurbFix;
}

void
BCInletTotalMassFlow::LocalFunc(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& U, const Block& block
    )
{
    double* Ui = U(iInterior);
    double* Ri = block.Radius()(iInterior);
    double* Rg = block.Radius()(iGhost);

    double gamma, gm1;
    gamma = Physics::GetInstance()->Gamma();
    gm1 = gamma - 1.0;

    // ASSUMPTION: Sn is normal to the entrainment velocity vector.

    // We assume that the velocity component normal to the boundary face stays the same across the boundary face.
    // So we know the normal component of the velocity in the ghost cell. We also know the direction of the velocity.
    // From these we can calcualte the ghost cell velocity, and along with total P and T we compute the rest.
    Vector3 ve(0.0, 0.0, 0.0); // entrainment velocity vector
    if (!block.IsStationary())
    {
        assert(block.IsRotating());
        RotationalMotion* rm = dynamic_cast<RotationalMotion*>(block.GetRigidBodyMotion());
        assert(rm != NULL);
        Vector3 pC(Ri[0], Ri[1], Ri[2]);
        ve = rm->EntrainmentVelocityAt(pC);
    }

    // Interior cell
    Vector3 n = Sn.Normalized();
    Vector3 vti;
    double vxi, vyi, vzi, rhoei, pi, ci, mdoti, vni, vti, vniabs;
    double romegaSqi, rhokeri;

    vxi = Ui[1] / Ui[0];
    vyi = Ui[2] / Ui[0];
    vzi = Ui[3] / Ui[0];
    romegaSqi = Ri[3];
    rhokeri = 0.5 * Ui[0] * romegaSqi;
    rhoei = Ui[4] - 0.5 * Ui[0] * (vxi * vxi + vyi * vyi + vzi * vzi) + rhokeri;
    pi = gm1 * rhoei;
    ci = std::sqrt(gamma * pi / Ui[0]);
    mdoti = Ui[1] * Sn.X() + Ui[2] * Sn.Y() + Ui[3] * Sn.Z();
    mMassFlowLocal += rhovni;

    vniabs = (vxi + ve.X()) * n[0] + (vyi + ve.Y()) * n[1] + (vzi + ve.Z()) * n[2];
    vni    = vxi * n[0] + vyi * n[1] + vzi * n[2];
    assert(vni >= 0.0);
    vti   = Vector3(vxi, vyi, vzi) - vni * n;

    // Compute what the static pressure and temperature should be based on the current mass flow rate w.r.t. the target.
    assert(mMassFlowFamilyCurrent > 0.0);
    const double SCALE_LOWER = 0.9, SCALE_UPPER = 1.1;
    double normalVelocityScale = std::max(SCALE_LOWER, std::min(SCALE_UPPER, mMassFlowFamilyTarget / mMassFlowFamilyCurrent));

    // Ghost cell velocity
    double vng;
    Vector3 vg, vgabs;
    vng = normalVelocityScale * vni;
    vg = vng * n + vti;
    vgabs = vg + ve;

    // From total temp & pres, static temp & pres are computed from isentropic relations.
    double Tg, Pg, rhog, rhoeg, romegaSqg, rhokerg, rhoetg;
    Tg = mTotalTemperature - 0.5 * gm1 * vgabs.MagSq() / gamma;
    Pg = mTotalPressure / std::pow(mTotalTemperature / Tg, gamma / gm1);
    rhog = Pg / Tg;
    rhoeg = Pg / gm1;
    romegaSqg = Rg[3];
    rhokerg = 0.5 * rhog * romegaSqg;
    rhoetg = rhoeg + 0.5 * vg.MagSq() * rhog - rhokerg;

    double UG[5];
    UG[0] = rhog;
    UG[1] = rhog * vg[0];
    UG[2] = rhog * vg[1];
    UG[3] = rhog * vg[2];
    UG[4] = rhoetg;

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
BCInletTotalMassFlow::LocalFuncTurb(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& UT, const Structured<double>& U, const Block& block
    )
{
    IndexIJK ijk = iGhost;
    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* ut = UT(ijk);
        double* u = U(ijk);
        double rho = u[0];
        // FIXME: model-specific
        ut[0] = rho * mTurbFix[0];
        ut[1] = rho * mTurbFix[1];
        ijk = ijk - deltaInterior;
    }
}

void
BCInletTotalMassFlow::PrepareStart(const Structured<double>& U)
{
    mMassFlowLocal = 0.0;
}

