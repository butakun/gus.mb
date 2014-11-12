// $Id: BCExtrapolate.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_BC_EXTRAPOLATE_H__
#define INCLUDED_BC_EXTRAPOLATE_H__

#include "BCPlanarLocal.h"

class BCExtrapolate : public BCPlanarLocal
{
public:
    BCExtrapolate(const IndexRange& range, Direction direction) : BCPlanarLocal(range, direction) {}
    virtual ~BCExtrapolate() {}

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

#endif // INCLUDED_BC_EXTRAPOLATE_H__

