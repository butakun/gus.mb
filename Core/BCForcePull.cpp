// $Id: BCForcePull.cpp 306 2013-10-02 07:03:25Z kato $

#include "Communicator.h"
#include "BCForcePull.h"
#include "Physics.h"

BCForcePull::BCForcePull(const FlowModel& model, const IndexRange& meshRange, Direction direction)
:   BCPlanarLocal(meshRange, direction),
    mModel(model)
{
}

void
BCForcePull::LocalFunc(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& U, const Block& block
    )
{
    double* Ui = U(iInterior);

    //double gamma;
    //gamma = Physics::GetInstance()->Gamma();

    double* Ri = block.Radius()(iInterior);
    //double* Rg = block.Radius()(iGhost);

    double vesqI, rhoetGhost, rhoVMagGhost;
    vesqI     = Ri[3]; // entrainment velocity squared at interior cell
    rhoetGhost = Ui[4] + 0.5 * Ui[0] * vesqI;
    rhoVMagGhost = std::sqrt(Ui[1] * Ui[1] + Ui[2] * Ui[2] + Ui[3] * Ui[3]);

    Vector3 nv = -Sn;
    nv.Normalize();

    Vector3 rhoVGhost = rhoVMagGhost * nv;

    double UGg[5], UG[5];
    UGg[0] = Ui[0];
    UGg[1] = rhoVGhost.X();
    UGg[2] = rhoVGhost.Y();
    UGg[3] = rhoVGhost.Z();
    UGg[4] = rhoetGhost;
    mModel.FromGlobalToLocal(UG, UGg, block, iGhost);

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
BCForcePull::LocalFuncTurb(
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
