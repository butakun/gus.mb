// $Id: BCViscousWall.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_BC_VISCOUS_WALL_H__
#define INCLUDED_BC_VISCOUS_WALL_H__

#include "BCPlanarLocal.h"

class BlockLocalFunctor
{
public:
    BlockLocalFunctor(Block& block) : mBlock(block) {}
    void operator () (
        const IndexIJK& iFace,
        const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
        const Vector3& Sn, const Block& block
        )
    {
        Structured<double> TurMuK = mBlock.TurMuK();
        double* turmukI = TurMuK(iInterior);
        double* turmukG = TurMuK(iGhost);
        turmukG[0] = -turmukI[0];
        turmukG[1] = -turmukI[1];
    }

private:
    Block& mBlock;
};

class BCViscousWall : public BCPlanarLocal
{
public:
    BCViscousWall(const IndexRange& meshRange, Direction direction, double TWall = -1.0);
    virtual ~BCViscousWall() {}

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
};

#endif // INCLUDED_BC_VISCOUS_WALL_H__

