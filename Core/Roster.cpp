/*
    gus.mb, an open source flow solver.
    Copyright (C) 2016 Hiromasa Kato <hiromasa at gmail.com>

    This file is part of gus.mb.

    gus.mb is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    gus.mb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Roster.h"
#include "Block.h"
#include "Communicator.h"
#include "BC.h"
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

const BlockPatch&
Roster::GetBlockPatch(int blockPatchID) const
{
    for (BlockPatchSet::const_iterator i = mBlockPatchSet.begin();
        i != mBlockPatchSet.end(); ++i)
    {
        if (i->UniqueID() == blockPatchID)
            return *i;
    }
    assert(false);
}

const BlockPatches&
Roster::RegisterBlockPatchFamily(int familyID, const BlockPatches& bps)
{
    assert(mBlockPatchFamilies.find(familyID) == mBlockPatchFamilies.end());
    return mBlockPatchFamilies[familyID] = bps;
}

const BlockPatches&
Roster::GetBlockPatchFamily(int familyID) const
{
    return mBlockPatchFamilies.at(familyID);
}

const char* MODELNAME = "CompressibleConservativeTwoEqTurb"; // FIXME: temporary

void
Roster::RegisterBC(int blockID, const Model* model, BC* bc)
{
#if 0
    mModelBlockBCs[model->Name()][blockID].push_back(bc);
#else
    mModelBlockBCs[MODELNAME][blockID].push_back(bc);
#endif
}

const std::vector<BC*>&
Roster::GetBCs(int blockID, const Model* model) const
{
    return mModelBlockBCs.find(MODELNAME)->second.find(blockID)->second;
}

void
Roster::ApplyBCs(const Model* model)
{
#if 0
    assert(mModelBlockBCs.find(model->Name()) != mModelBlockBCs.end());
    BlockBCs blockBCs = mModelBlockBCs[model->Name()];
#else
    assert(mModelBlockBCs.find(MODELNAME) != mModelBlockBCs.end());
    BlockBCs blockBCs = mModelBlockBCs[MODELNAME];
#endif

    for (BlockBCs::iterator i = blockBCs.begin(); i != blockBCs.end(); ++i)
    {
        int blockID = i->first;
        BCs bcs = i->second;
        Block* block = dynamic_cast<Block*>(GetBlock(blockID));
        assert(block != NULL);
        for (BCs::iterator j = bcs.begin(); j != bcs.end(); ++j)
        {
            (*j)->Apply(*block, block->U());
        }
        block->FillCornerGhosts();
    }
}

void
Roster::ApplyTurbBCs(const Model* model)
{
#if 0
    assert(mModelBlockBCs.find(model->Name()) != mModelBlockBCs.end());
    BlockBCs blockBCs = mModelBlockBCs[model->Name()];
#else
    assert(mModelBlockBCs.find(MODELNAME) != mModelBlockBCs.end());
    BlockBCs blockBCs = mModelBlockBCs[MODELNAME];
#endif

    for (BlockBCs::iterator i = blockBCs.begin(); i != blockBCs.end(); ++i)
    {
        int blockID = i->first;
        BCs bcs = i->second;
        Block* block = dynamic_cast<Block*>(GetBlock(blockID));
        assert(block != NULL);
        for (BCs::iterator j = bcs.begin(); j != bcs.end(); ++j)
        {
            (*j)->ApplyTurb(*block, block->UT(), block->U());
        }
        block->FillCornerGhosts();
    }
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

#if 0
const Structured<double>&
Roster::GetMeshForBlockPatch(const BlockPatch& bp)
{
    if (mBlockPatchMeshes.find(bp.UniqueID()) == mBlockPatchMeshes.end())
    {
        const Structured<double>& XYZ = GetMeshForBlock(bp.BlockID());
        assert(false); // FIXME
    }
}
#endif

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

