// $Id: BCExtrapolate.cpp 277 2013-06-04 01:58:51Z kato $

#include "BCExtrapolate.h"

void
BCExtrapolate::LocalFunc(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& U, const Block& block
    )
{
    IndexIJK iI = iInterior;
    IndexIJK iG = iGhost;
    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* UI = U(iI);
        double* UG = U(iG);
        for (int l = 0; l < U.DOF(); ++l)
        {
            UG[l] = UI[l];
        }
        iG -= deltaInterior;
    }
}

void
BCExtrapolate::LocalFuncTurb(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& UT, const Structured<double>& U, const Block& block
    )
{
    IndexIJK iI = iInterior;
    IndexIJK iG = iGhost;
    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* UI = UT(iI);
        double* UG = UT(iG);
        for (int l = 0; l < UT.DOF(); ++l)
        {
            UG[l] = UI[l];
        }
        iG -= deltaInterior;
    }
}

