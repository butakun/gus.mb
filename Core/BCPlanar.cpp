// $Id: BCPlanar.cpp 231 2012-05-11 01:11:28Z kato $

#include "BCPlanar.h"
#include "IndexUtils.h"
#include <cassert>

IndexRange
BCPlanar::RindRange() const
{
    return IndexUtils::FromMeshRangeToCellRange(MeshRange(), PatchDirection(), 1); // FIXME: only single ghost layer?
}

void
BCPlanar::SetMaskWithAValue(Structured<int>& mask, int value) const
{
    IndexRange rindRange = RindRange();
    for (IndexIterator itor(rindRange); !itor.IsEnd(); itor.Advance())
    {
        IndexIJK ijk = itor.Index();
        *mask(ijk) = value;
    }
}

