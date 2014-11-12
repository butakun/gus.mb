// $Id: BCSymmetry.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_BC_SYMMETRY_H__
#define INCLUDED_BC_SYMMETRY_H__

#include "BCPlanarLocal.h"

class BCSymmetry : public BCPlanarLocal
{
public:
    BCSymmetry(const IndexRange& meshRange, Direction direction)
    : BCPlanarLocal(meshRange, direction) {}

    virtual ~BCSymmetry() {}

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
};

#endif // INCLUDED_BC_SYMMETRY_H__

