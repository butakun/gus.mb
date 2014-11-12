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

    //double* Ri = block.Radius()(iInterior);
    double* Rg = block.Radius()(iGhost);

    double rhoGhost, uGhost, vGhost, wGhost, pGhost, rhoke, rhoetGhost;
    double romegaSq, rhoker;
    rhoGhost = Ui[0];
    uGhost = Ui[1] / Ui[0];
    vGhost = Ui[2] / Ui[0];
    wGhost = Ui[3] / Ui[0];
    pGhost = mStaticPressure;
    romegaSq = Rg[3];
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
