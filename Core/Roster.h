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
#ifndef INCLUDE_ROSTER_H__
#define INCLUDE_ROSTER_H__

#include "Communicator.h"
#include "VirtualBlock.h"
#include "Block.h"
#include "Structured.h"
#include "AbuttingInterface.h"
#include <set>
#include <map>
#include <vector>

class BC;
class Model;

class Roster
{
public:
    typedef std::set<BlockPatch> BlockPatchSet;
    typedef std::map<int, BlockPatches> BlockPatchesMap;

    static Roster* GetInstance();

    virtual ~Roster();

    void RegisterBlock(int rank, VirtualBlock* block);
    int GetRankOf(int blockID) const;
    VirtualBlock* GetBlock(int blockID) const;
    int GetNumberOfBlocks() const { return mBlockMap.size(); }

    // BlockPatch database
    const BlockPatch& RegisterBlockPatch(const BlockPatch& bp);
    const BlockPatch& GetBlockPatch(int blockID, const IndexRange& meshRange) const;
    const BlockPatch& GetBlockPatch(int blockPatchID) const;
    const BlockPatches& RegisterBlockPatchFamily(int familyID, const BlockPatches& bps);
    const BlockPatches& GetBlockPatchFamily(int familyID) const;
    const BlockPatchesMap& BlockPatchFamilies() const { return mBlockPatchFamilies; }

    // BCs
    void RegisterBC(int blockID, const Model* model, BC* bc);
    const std::vector<BC*>& GetBCs(int blockID, const Model* model) const;
    void ApplyBCs(const Model* model);
    void ApplyTurbBCs(const Model* model);

    // Interfaces
    void AddInterface(AbuttingInterface* interface);
    const Interfaces& GetInterfaces() const { return mInterfaces; }
    Interfaces& GetInterfaces() { return mInterfaces; }

    // Global mesh data (temporary)
    void SetMeshFileName(const char* filename, double scaling);
    const Structured<double>& GetMeshForBlock(int blockID);
    const Structured<double>& GetMeshForBlockPatch(const BlockPatch& bp);
    void ReleaseCachedMesh();

    // a tag is a positive (non-zero) integer.
    int CreateNewTag(const char* name);
    int FindTagByName(const char* name) const; // returns zero if the tag does not exist.
    template <class T> void AddBlockData(int blockID, int tag, const char* name, Structured<T>& Data);
    template <class T> const Structured<T>& GetData(int blockID, int tag) const;

    Structured<double>& GetBlockData(Block& block, const char* name) const; // FIXME
    Structured<double>& GetBlockData(int blockID, const char* name) const; // FIXME

protected:
    Roster();

private:
    static Roster* Instance;

    // Blocks
    std::map<int, int> mBlockRankMap;
    std::map<int, VirtualBlock*> mBlockMap;

    // BlockPatches
    BlockPatchSet mBlockPatchSet;
    BlockPatchesMap mBlockPatchFamilies;

    // Interfaces
    Interfaces mInterfaces;

    // Cached mesh (separate from the mesh retained by local blocks, used for interfaces, non-matching patches, etc.)
    std::string mMeshFileName;
    double mMeshScaling;
    std::map<int, Structured<double> > mCachedMesh;
    std::map<int, Structured<double> > mBlockPatchMeshes;

    // BCs
    typedef std::vector<BC*> BCs;
    typedef std::map<int, BCs> BlockBCs;
    typedef std::map<std::string, BlockBCs> ModelBlockBCs;
    ModelBlockBCs mModelBlockBCs;

    // BC Family data
    class BCFamilyData
    {
    public:
        int BlockPatchID;
        double Data;
        BCFamilyData() : BlockPatchID(0), Data(0.0) {}
    };

    typedef std::map<int, std::map<std::string, BCFamilyData> > BCFamilyDataMap;

    // Tagged BlockData
    class TaggedBlockDataBase
    {
    public:
        int BlockID;
        int Tag;
        std::string Name;
        TaggedBlockDataBase()
        : BlockID(-1), Tag(-1)
        {
        }
        TaggedBlockDataBase(int blockID, int tag, const char* name)
        : BlockID(blockID), Tag(tag), Name(name)
        {
        }
    };

    template <class T>
    class TaggedBlockData : public TaggedBlockDataBase
    {
    public:
        Structured<T> Data;
        TaggedBlockData()
        {
        }
        TaggedBlockData(int blockID, int tag, const char* name)
        : TaggedBlockDataBase(blockID, tag, name)
        {
        }
        TaggedBlockData(int blockID, int tag, const char* name, Structured<T>& data)
        : TaggedBlockDataBase(blockID, tag, name), Data(data)
        {
        }
        bool operator == (const TaggedBlockData<T>& td) const
        {
            return BlockID == td.BlockID && Tag == td.Tag;
        }
    };

    typedef std::map<std::string, int> NameTagMap;
    NameTagMap mNameTagMap;

    typedef std::map<int, TaggedBlockDataBase*> TaggedBlockDataMap;
    typedef std::map<int, TaggedBlockDataMap> BlockTagMap;
    BlockTagMap mBlockTagMap;
};

template <class T>
void
Roster::AddBlockData(int blockID, int tag, const char* name, Structured<T>& data)
{
    Communicator::GetInstance()->Console() << "Roster::AddBlockData: block " << blockID << " tag = " << tag << " name = " << name << std::endl;
    TaggedBlockDataMap& tbdm = mBlockTagMap[blockID];
    TaggedBlockDataMap::const_iterator i = tbdm.find(tag);
    if (i != tbdm.end())
    {
        // Data with the same signature exists.
        Communicator::GetInstance()->Console() << "Roster::AddBlockData:tag already exists." << std::endl;
        throw 666;
    }
    tbdm[tag] = new TaggedBlockData<T>(blockID, tag, name, data);
}

template <class T>
const Structured<T>&
Roster::GetData(int blockID, int tag) const
{
    BlockTagMap::const_iterator i = mBlockTagMap.find(blockID);
    if (i == mBlockTagMap.end())
    {
        Communicator::GetInstance()->Console() << "Roster::GetData:couldn't find data for block ID " << blockID << std::endl;
        throw 666;
    }

    const TaggedBlockDataMap& tdm = i->second;
    TaggedBlockDataMap::const_iterator j = tdm.find(tag);
    if (j == tdm.end())
    {
        Communicator::GetInstance()->Console() << "Roster::GetData:couldn't find data with tag = " << tag << " in block " << blockID << std::endl;
        throw 666;
    }

    TaggedBlockData<T>* tbd = dynamic_cast<TaggedBlockData<T>*>(j->second);
    assert(tbd != NULL);

    return tbd->Data;
}

#endif // INCLUDE_ROSTER_H__

