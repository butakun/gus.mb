// $Id: BCCopy.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_BC_COPY_H__
#define INCLUDED_BC_COPY_H__

#include "BCPlanarLocal.h"

class BCCopy : public BCPlanarLocal
{
public:
    BCCopy(const IndexRange& range, Direction direction) : BCPlanarLocal(range, direction) {}
    virtual ~BCCopy() {}

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

#endif // INCLUDED_BC_COPY_H__

