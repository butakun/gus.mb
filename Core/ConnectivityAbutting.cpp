// $Id: ConnectivityAbutting.cpp 41 2010-07-06 10:48:07Z kato $

#include "ConnectivityAbutting.h"

ConnectivityAbutting::ConnectivityAbutting(
    const IndexRange& meshRange, Direction dir, const Block& block,
    const IndexRange& donorMeshRange, int donorBlockID,
    int tag
    )
:   Connectivity(meshRange, dir),
    mBlockID(block.ID()),
    mDonorBlockID(donorBlockID),
    mTag(tag)
{
}

ConnectivityAbutting::~ConnectivityAbutting()
{
}

