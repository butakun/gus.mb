// $Id: BCSymmetry.cpp 277 2013-06-04 01:58:51Z kato $

#include "BCSymmetry.h"
#include "Vector3.h"

void
BCSymmetry::LocalFunc(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& U, const Block& block
    )
{
    IndexIJK iI = iInterior;
    IndexIJK iG = iGhost;
    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* UGhost = U(iG);
        double* UInterior = U(iI);

        for (int l = 0; l < U.DOF(); ++l)
        {
            UGhost[l] = UInterior[l];
        }

        Vector3 V(UInterior[1], UInterior[2], UInterior[3]);
        V = V / UInterior[0];

        Vector3 n = Sn / Sn.Mag();
        Vector3 Vn = n * dot_product(V, n);
        Vector3 Vghost = V - 2.0 * Vn;
        UGhost[1] = UGhost[0] * Vghost[0];
        UGhost[2] = UGhost[0] * Vghost[1];
        UGhost[3] = UGhost[0] * Vghost[2];

        iI += deltaInterior;
        iG -= deltaInterior;
    }
}

void
BCSymmetry::LocalFuncTurb(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& UT, const Structured<double>& U, const Block& block
    )
{
    IndexIJK iI = iInterior;
    IndexIJK iG = iGhost;
    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* UGhost = UT(iG);
        double* UInterior = UT(iI);

        for (int l = 0; l < UT.DOF(); ++l)
        {
            UGhost[l] = UInterior[l];
        }

        iI += deltaInterior;
        iG -= deltaInterior;
    }
}

