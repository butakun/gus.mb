// $Id: BCFix.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_BC_FIX_H__
#define INCLUDED_BC_FIX_H__

#include "BCPlanarLocal.h"

class TurbulenceSpec;

class BCFix : public BCPlanarLocal
{
public:
    BCFix(const IndexRange& range, Direction direction, int ndof, double* ufix, int ndoft, const TurbulenceSpec* turbSpec);
    virtual ~BCFix();

    double* U() const { return UFix; }
    void SetTo(int ndof, double* u, int ndoft, const TurbulenceSpec* turbSpec);

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
    double* UFix;
    double* TurbFix;
};

#endif // INCLUDED_BC_FIX_H__

