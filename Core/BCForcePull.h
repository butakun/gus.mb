// $Id: BCForcePull.h 291 2013-07-25 10:25:11Z kato $
#ifndef INCLUDE_BC_FORCE_PULL_H__
#define INCLUDE_BC_FORCE_PULL_H__

#include "BCPlanarLocal.h"
#include "FlowModel.h"

class BCForcePull : public BCPlanarLocal
{
public:
    enum Axis { X, Y, Z };

    // The velocity vector is forced to be along the unit vector (vaxial, vradial).
    BCForcePull(const FlowModel& model, const IndexRange& meshRange, Direction direction);
    virtual ~BCForcePull() {}

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
    const FlowModel& mModel;
};

#endif // INCLUDE_BC_FORCE_PULL_H__

