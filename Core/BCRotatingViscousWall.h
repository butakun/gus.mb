// $Id: BCViscousWall.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_BC_ROTATING_VISCOUS_WALL_H__
#define INCLUDED_BC_ROTATING_VISCOUS_WALL_H__

#include "BCPlanarLocal.h"
#include "BCViscousWall.h"
#include "RigidBodyMotion.h"

class BCRotatingViscousWall : public BCPlanarLocal
{
public:
    BCRotatingViscousWall(const IndexRange& meshRange, Direction direction, double TWall, const Vector3& angularVelocity);
    virtual ~BCRotatingViscousWall();

    virtual void SetMask(Structured<int>& mask) const;

    void ApplyTurMuK(Block& block)
    {
        BlockLocalFunctor ftor(block);
        Apply<BlockLocalFunctor>(block, ftor);
    }

protected:
    virtual void LocalFunc(
        const IndexIJK& iFace,
        const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
        const Vector3& Sn, Structured<double>& U, const Block& block
        );
    virtual void LocalFuncTurb(
        const IndexIJK& iFace,
        const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
        const Vector3& Sn, Structured<double>& UT, const Structured<double>& U, const Block& block
        );

private:
    double mTWall;
    RotationalMotion* mRBM;
};

#endif // INCLUDED_BC_ROTATING_VISCOUS_WALL_H__

