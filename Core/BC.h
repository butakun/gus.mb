// $Id: BC.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_BC_H__
#define INCLUDED_BC_H__

#include "Block.h"
#include "Structured.h"

class BC
{
public:
    BC(const IndexRange& meshRange);
    virtual ~BC() {}

    virtual void Apply(const Block& block, Structured<double>& U) = 0;
    virtual void ApplyTurb(const Block& block, Structured<double>& UT, const Structured<double>& U) = 0;
    virtual void SetMask(Structured<int>& mask) const {}

    const IndexRange& MeshRange() const { return mMeshRange; }

protected:

private:
    IndexRange mMeshRange;
};

#endif // INCLUDED_BC_H__

