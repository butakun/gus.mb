// $Id: BCInletTotal.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDE_BC_INLET_TOTAL_H__
#define INCLUDE_BC_INLET_TOTAL_H__

#include "BCPlanarLocal.h"

class TurbulenceSpec;

class BCInletTotal : public BCPlanarLocal
{
public:
    BCInletTotal(
        const IndexRange& meshRange, Direction direction,
        double pt, double tt, const Vector3& vdir, int ndoft, const TurbulenceSpec* turbSpec
        );
    virtual ~BCInletTotal();

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
    double mTotalPressure;
    double mTotalTemperature;
    Vector3 mVdir;
    double* mTurbFix;
};

#endif // INCLUDE_BC_INLET_TOTAL_H__

