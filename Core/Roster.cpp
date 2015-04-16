// $Id: Roster.cpp 306 2013-10-02 07:03:25Z kato $

#include "Roster.h"
#include "Block.h"
#include "Communicator.h"
#include "CGNSReader.h"
#include <algorithm>
#include <cassert>

Roster* Roster::Instance = NULL;

Roster* Roster::GetInstance()
{
    if (Instance == NULL)
    {
        Instance = new Roster();
    }
    return Instance;
}

Roster::Roster()
{
}

Roster::~Roster()
{
}

void
Roster::RegisterBlock(int rank, VirtualBlock* block)
{
    mBlockRankMap[block->ID()] = rank;
    mBlockMap[block->ID()] = block;
}

int
Roster::GetRankOf(int blockID) const
{
    std::map<int, int>::const_iterator i = mBlockRankMap.find(blockID);
    assert(i != mBlockRankMap.end());
    return i->second;
}

VirtualBlock*
Roster::GetBlock(int blockID) const
{
    std::map<int, VirtualBlock*>::const_iterator i = mBlockMap.find(blockID);
    assert(i != mBlockMap.end());
    return i->second;
}

const BlockPatch&
Roster::RegisterBlockPatch(const BlockPatch& bp)
{
    assert(mBlockPatchSet.find(bp) == mBlockPatchSet.end());
    BlockPatchSet::const_iterator i = mBlockPatchSet.insert(bp).first;
    const BlockPatch& bp2 = *i;
    return bp2;
}

const BlockPatch&
Roster::GetBlockPatch(int blockID, const IndexRange& meshRange) const
{
    BlockPatch bp(blockID, meshRange);
    BlockPatchSet::const_iterator i = mBlockPatchSet.find(bp);
    assert(i != mBlockPatchSet.end());
    return *i;
}

void
Roster::AddInterface(AbuttingInterface* interface)
{
    mInterfaces.push_back(interface);
}

void
Roster::SetMeshFileName(const char* filename, double scaling)
{
    // FIXME: mutex!
    mMeshFileName = filename;
    mMeshScaling = scaling;
    ReleaseCachedMesh();
}

const Structured<double>&
Roster::GetMeshForBlock(int blockID)
{
    if (mCachedMesh.find(blockID) == mCachedMesh.end())
    {
        CGNSReader* reader = new CGNSReader(mMeshFileName.c_str());
        Structured<double> XYZ;
        std::string zoneName;
        int Z = blockID;
        reader->ReadMesh(Z, XYZ, zoneName);
        delete reader;

        mCachedMesh[blockID] = XYZ;
    }

    return mCachedMesh[blockID];
}

void
Roster::ReleaseCachedMesh()
{
    for (std::map<int, Structured<double> >::const_iterator i = mCachedMesh.begin();
        i != mCachedMesh.end(); ++i)
    {
        delete i->second.Data;
    }
    mCachedMesh.clear();
}

int
Roster::CreateNewTag(const char* name)
{
    int tag;

    std::string tagName(name);
    NameTagMap::const_iterator i = mNameTagMap.find(tagName);
    if (i == mNameTagMap.end())
    {
        // the tag with this name doesn't exist, so create a new one.
        tag = mNameTagMap.size() + 1;
        mNameTagMap[name] = tag;
        Communicator::GetInstance()->Console() << "Roster::CreateNewTag: new tag created " << tagName << std::endl;
    }
    else
    {
        tag = i->second;
    }

    return tag;
}

int
Roster::FindTagByName(const char* name) const
{
    std::string tagName(name);
    NameTagMap::const_iterator i = mNameTagMap.find(tagName);
    if (i == mNameTagMap.end())
    {
        Communicator::GetInstance()->Console() << "Roster::FindTagByName: no tag named " << tagName << std::endl;
        return 0;
    }

    return i->second;
}

Structured<double>&
Roster::GetBlockData(Block& block, const char* name) const
{
    std::string nam(name);
    if (nam == "U")
        return block.U();
    else if (nam == "UT")
        return block.UT();
    assert(false);
    return block.U(); // FIXME: gcc complains unless I return something.
}

Structured<double>&
Roster::GetBlockData(int blockID, const char* name) const
{
    VirtualBlock* vb = GetBlock(blockID);
    Block* block = dynamic_cast<Block*>(vb);
    assert(block != NULL);

    return GetBlockData(*block, name);
}

