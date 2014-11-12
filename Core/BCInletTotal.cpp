// $Id: BCInletTotal.cpp 306 2013-10-02 07:03:25Z kato $

#include "Communicator.h"
#include "BCInletTotal.h"
#include "Physics.h"
#include "TurbulenceSpec.h"
#include "RigidBodyMotion.h"
#include <cassert>

BCInletTotal::BCInletTotal(
    const IndexRange& meshRange, Direction direction,
    double pt, double tt, const Vector3& vdir,
    int ndoft, const TurbulenceSpec* turbSpec
    )
:   BCPlanarLocal(meshRange, direction),
    mTotalPressure(pt / Physics::GetInstance()->PRef()),
    mTotalTemperature(tt / Physics::GetInstance()->TRef()),
    mVdir(vdir),
    mTurbFix(new double[ndoft])
{
    mVdir.Normalize();

    assert(ndoft == 2);

    // FIXME: how do we make this model-neutral?
    mTurbFix[0] = turbSpec->TKE_Nondimensional();
    mTurbFix[1] = turbSpec->Omega_Nondimensional();

    Communicator::GetInstance()->Console() << "BCInletTotal: pt(nondim) = " << mTotalPressure << ", Tt(nondim) = " << mTotalTemperature << std::endl;
}

BCInletTotal::~BCInletTotal()
{
    delete[] mTurbFix;
}

void
BCInletTotal::LocalFunc(
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

    Vector3 n = Sn.Normalized();
    double vxi, vyi, vzi, vni, rhoei, pi, ci;
    double romegaSqi, rhokeri;

    // vxi, vyi, vzi are relative velocity components, whereas vni is in absolute frame.
    vxi = Ui[1] / Ui[0];
    vyi = Ui[2] / Ui[0];
    vzi = Ui[3] / Ui[0];
    romegaSqi = Ri[3];
    rhokeri = 0.5 * Ui[0] * romegaSqi;
    rhoei = Ui[4] - 0.5 * Ui[0] * (vxi * vxi + vyi * vyi + vzi * vzi) + rhokeri;
    pi = gm1 * rhoei;
    ci = std::sqrt(gamma * pi / Ui[0]);

    // vni is in absolute frame.
    vni = (vxi + ve.X()) * n[0] + (vyi + ve.Y()) * n[1] + (vzi + ve.Z()) * n[2];
    //assert(vni >= 0.0);
    vni = std::abs(vni);

    double ndotVdir = dot_product(n, mVdir); // FIXME: shouldn't this be written as Vdir dot n? though it's the same value.
    double Vabs = vni / ndotVdir;

    // Vrel = ghost cell velocity in relative frame
    Vector3 Vrel = Vabs * mVdir - ve; // FIXME: assumes the entrainment velocity in the ghost cell is the same as the interior cell. wouldn't work for centrifugal turbine inlet.

    double Mg = Vabs / ci;
    double isen = 1.0 + 0.5 * gm1 * Mg * Mg;
    double Pg, Tg, rhog, rhoeg, rhoetg;
    double romegaSqg, rhokerg;
    Tg = mTotalTemperature / isen;
    Pg = mTotalPressure / std::pow(isen, gamma / gm1);
    rhog = Pg / Tg;
    rhoeg = Pg / gm1;
    romegaSqg = Rg[3];
    rhokerg = 0.5 * rhog * romegaSqg;

    rhoetg = rhoeg + 0.5 * Vrel.MagSq() * rhog - rhokerg;

    Vector3 Vg = Vrel;

    double UG[5];
    UG[0] = rhog;
    UG[1] = rhog * Vg[0];
    UG[2] = rhog * Vg[1];
    UG[3] = rhog * Vg[2];
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
BCInletTotal::LocalFuncTurb(
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
