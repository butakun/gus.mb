// $Id: IndexUtils.h 207 2012-03-07 03:42:26Z kato $
#ifndef INCLUDED_INDEX_UTILS_H__
#define INCLUDED_INDEX_UTILS_H__

#include "Structured.h"

enum Direction { I = 1, J = 2, K = 3, INEG = -1, JNEG = -2, KNEG = -3 };

namespace IndexUtils
{

inline
Direction Opposite(Direction dir)
{
    return Direction(-int(dir));
}

inline
IndexRange FromMeshRangeToMetricRange(const IndexRange& meshRange, Direction direction)
{
    assert(meshRange.IsCanonical());
    switch (direction)
    {
    case I:
    case INEG:
        assert(meshRange.Start.I == meshRange.End.I);
        return IndexRange(
            meshRange.Start.I, meshRange.Start.J + 1, meshRange.Start.K + 1,
            meshRange.End.I, meshRange.End.J, meshRange.End.K
            );
    case J:
    case JNEG:
        assert(meshRange.Start.J == meshRange.End.J);
        return IndexRange(
            meshRange.Start.I + 1, meshRange.Start.J, meshRange.Start.K + 1,
            meshRange.End.I, meshRange.End.J, meshRange.End.K
            );
    case K:
    case KNEG:
        assert(meshRange.Start.K == meshRange.End.K);
        return IndexRange(
            meshRange.Start.I + 1, meshRange.Start.J + 1, meshRange.Start.K,
            meshRange.End.I, meshRange.End.J, meshRange.End.K
            );
    }

    // you should not be here
    assert(false);
    return meshRange;
}

inline
IndexRange FromMeshRangeToCellRange(const IndexRange& meshRange, Direction direction, int nGhostLayers)
{
    assert(meshRange.IsCanonical());
    switch (direction)
    {
    case I:
        assert(meshRange.Start.I == meshRange.End.I);
        return IndexRange(
            meshRange.Start.I + 1, meshRange.Start.J + 1, meshRange.Start.K + 1,
            meshRange.End.I + nGhostLayers, meshRange.End.J, meshRange.End.K
            );
    case INEG:
        assert(meshRange.Start.I == meshRange.End.I);
        return IndexRange(
            meshRange.Start.I - (nGhostLayers - 1), meshRange.Start.J + 1, meshRange.Start.K + 1,
            meshRange.End.I, meshRange.End.J, meshRange.End.K
            );
    case J:
        assert(meshRange.Start.J == meshRange.End.J);
        return IndexRange(
            meshRange.Start.I + 1, meshRange.Start.J + 1, meshRange.Start.K + 1,
            meshRange.End.I, meshRange.End.J + nGhostLayers, meshRange.End.K
            );
    case JNEG:
        assert(meshRange.Start.J == meshRange.End.J);
        return IndexRange(
            meshRange.Start.I + 1, meshRange.Start.J - (nGhostLayers - 1), meshRange.Start.K + 1,
            meshRange.End.I, meshRange.End.J, meshRange.End.K
            );
    case K:
        assert(meshRange.Start.K == meshRange.End.K);
        return IndexRange(
            meshRange.Start.I + 1, meshRange.Start.J + 1, meshRange.Start.K + 1,
            meshRange.End.I, meshRange.End.J, meshRange.End.K + nGhostLayers
            );
    case KNEG:
        assert(meshRange.Start.K == meshRange.End.K);
        return IndexRange(
            meshRange.Start.I + 1, meshRange.Start.J + 1, meshRange.Start.K - (nGhostLayers - 1),
            meshRange.End.I, meshRange.End.J, meshRange.End.K
            );
    }

    // you should not be here
    assert(false);
    return meshRange;
}

Direction PatchDirection(const IndexRange& blockMeshRange, const IndexRange& patchMeshRange);

Direction PatchDirection(const IndexRange& patchMeshRange); // only returns I, J, or K, not INEG, JNEG, KNEG.

}

#endif // INCLUDED_INDEX_UTILS_H__

