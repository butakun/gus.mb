// $Id: BCOutletStaticPressure.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDE_BC_OUTLET_STATIC_PRESSURE_H__
#define INCLUDE_BC_OUTLET_STATIC_PRESSURE_H__

#include "BCPlanarLocal.h"

class BCOutletStaticPressure : public BCPlanarLocal
{
public:
    BCOutletStaticPressure(const IndexRange& meshRange, Direction direction, double pressure);
    virtual ~BCOutletStaticPressure() {}

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
    double mStaticPressure;
};

#endif // INCLUDE_BC_OUTLET_STATIC_PRESSURE_H__

