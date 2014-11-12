// $Id: BCPlanar.h 231 2012-05-11 01:11:28Z kato $
#ifndef INCLUDED_BC_PLANAR_H__
#define INCLUDED_BC_PLANAR_H__

#include "BC.h"
#include "IndexUtils.h"

class BCPlanar : public BC
{
public:
    BCPlanar(const IndexRange& meshRange, Direction direction)
    : BC(meshRange), mDirection(direction) {}
    virtual ~BCPlanar() {}

    Direction PatchDirection() const { return mDirection; }
    IndexRange MetricRange() const { return IndexUtils::FromMeshRangeToMetricRange(MeshRange(), PatchDirection()); }
    IndexRange RindRange() const;
    Structured<double> Surface(const Block& block) const
    {
        switch (PatchDirection())
        {
        case I:
        case INEG:
            return block.Sxi();
        case J:
        case JNEG:
            return block.Seta();
        case K:
        case KNEG:
            return block.Szeta();
        }
        assert(false);
        return block.Sxi();
    }

protected:
    void SetMaskWithAValue(Structured<int>& mask, int value) const;

private:
    Direction mDirection;
};

#endif // INCLUDED_BC_PLANAR_H__

