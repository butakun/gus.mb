// $Id: BCInviscidWall.h 32 2010-05-28 10:37:03Z kato $
#ifndef INCLUDED_BC_INVISCID_WALL_H__
#define INCLUDED_BC_INVISCID_WALL_H__

#include "BCPlanarLocal.h"

#if 0
class BCInviscidWall : public BCPlanarLocal
{
public:
    BCInviscidWall(const IndexRange& meshRange, Direction direction);
    virtual ~BCInviscidWall();

    virtual void Apply(const Block& block, Structured<double>& U);

protected:

private:
};
#else

#include "Vector3.h"

inline
void BCInviscidWallFunc(
    const IndexIJK& iFace, const IndexIJK& dGhost, const IndexIJK& dInterior,
    const Vector3& Sn, Structured<double>& U, const Block& block
    )
{
    IndexIJK iGhost = iFace + dGhost;
    IndexIJK iInterior = iFace + dInterior;

    double* UGhost = U(iGhost);
    double* UInterior = U(iInterior);

    for (int l = 0; l < 5; ++l)
    {
        UGhost[l] = UInterior[l];
    }

    Vector3 V(UInterior[1], UInterior[2], UInterior[3]);
    V = V / UInterior[0];

    Vector3 n = Sn / Sn.Mag();
    Vector3 Vn = n * dot_product(V, n);
    Vector3 Vghost = V - 2.0 * Vn;
    UGhost[1] = UGhost[0] * Vghost[0];
    UGhost[2] = UGhost[1] * Vghost[1];
    UGhost[3] = UGhost[2] * Vghost[2];
}

typedef BCPlanarLocal<BCInviscidWallFunc> BCInviscidWall;
#endif

#endif // INCLUDED_BC_INVISCID_WALL_H__

