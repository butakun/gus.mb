// $Id: IndexUtils.cpp 166 2011-12-13 04:52:50Z kato $

#include "IndexUtils.h"

Direction
IndexUtils::PatchDirection(const IndexRange& blockMeshRange,  const IndexRange& patchMeshRange)
{
    IndexIJK shape = patchMeshRange.Shape();
    if (shape.I == 1)
    {
        if (patchMeshRange.Start.I == blockMeshRange.Start.I)
            return I;
        else if (patchMeshRange.End.I == blockMeshRange.End.I)
            return INEG;
        else
            assert(false);
    }
    else if (shape.J == 1)
    {
        if (patchMeshRange.Start.J == blockMeshRange.Start.J)
            return J;
        else if (patchMeshRange.End.J == blockMeshRange.End.J)
            return JNEG;
        else
            assert(false);
    }
    else if (shape.K == 1)
    {
        if (patchMeshRange.Start.K == blockMeshRange.Start.K)
            return K;
        else if (patchMeshRange.End.K == blockMeshRange.End.K)
            return KNEG;
        else
            assert(false);
    }
    else
    {
        std::cerr << "Patch shape is not planar: " << patchMeshRange << ", shape = " << shape << std::endl;
        assert(false);
    }

    // Should not be here
    throw 666;
}

Direction
IndexUtils::PatchDirection(const IndexRange& patchMeshRange)
{
    IndexIJK shape = patchMeshRange.Shape();
    if (shape.I == 1)
    {
        return I;
    }
    else if (shape.J == 1)
    {
        return J;
    }
    else if (shape.K == 1)
    {
        return K;
    }

    // Should not be here
    throw 666;
}

