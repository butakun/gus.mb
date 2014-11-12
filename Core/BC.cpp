// $Id: BC.cpp 30 2010-05-24 00:47:45Z kato $

#include "BC.h"
#include <cassert>

BC::BC(const IndexRange& meshRange)
:   mMeshRange(meshRange)
{
    assert(mMeshRange.IsCanonical());
}

