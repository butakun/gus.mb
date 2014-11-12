// $Id: BlockPatch.h 254 2012-10-07 17:51:59Z kato $
#ifndef INCLUDED_BLOCK_PATCH_H__
#define INCLUDED_BLOCK_PATCH_H__

#include "Structured.h"
#include "IndexUtils.h"
#include <iostream>

class BlockPatch
{
public:
    static BlockPatch New(int blockID, const IndexRange& meshRange, int uniqueID);

    int BlockID() const { return mBlockID; }
    const IndexRange& MeshRange() const { return mMeshRange; }
    const IndexRange& CellRange() const { return mCellRange; }
    IndexRange CellFaceRange() const;

    // returns the cell range of the first rind layer. the second (deeper) rind layer can be had by incrementing the rage by di3.
    void GhostCellRange(IndexRange& gcr, IndexIJK& i1, IndexIJK& i2, IndexIJK& di3) const;

    const IndexIJK& I1() const { return mI1; }
    const IndexIJK& I2() const { return mI2; }
    const IndexIJK& I3Ghost() const { return mI3Ghost; }
    const IndexIJK& I3Interior() const { return mI3Interior; }
    const IndexIJK& DI3Ghost() const { return mDI3Ghost; }

    int UniqueID() const
    {
#if 0
        int id = 0;
        IndexRange mr = MeshRange().Canonical();
        IndexIJK s = mr.Start;
        int d; // index normal to the patch surface. I -> 1, J -> 2, K -> 3
        IndexIJK shape = mr.Shape();
        if (shape.I == 1)
            d = 1;
        else if (shape.J == 1)
            d = 2;
        else if (shape.K == 1)
            d = 3;
        else
            throw 666;
        id = 1e8 * BlockID() + 1e7 * d + 1e6 * (s.K % 1000) + 1e3 * (s.J % 1000) + s.I % 1000;
        return id;
#else
        return mUniqueID;
#endif
    }

    bool operator == (const BlockPatch& bp) const
    {
        return BlockID() == bp.BlockID() && MeshRange().Canonical() == bp.MeshRange().Canonical();
    }

    // This needs to be implemented for this class to be used as a key to std::map.
    bool operator < (const BlockPatch& bp) const
    {
        if (BlockID() < bp.BlockID())
            return true;
        else if (BlockID() > bp.BlockID())
            return false;

        IndexRange r1 = MeshRange().Canonical();
        IndexRange r2 = bp.MeshRange().Canonical();
        if (r1.Start.I < r2.Start.I) return true;
        else if (r1.Start.I > r2.Start.I) return false;
        else
            if (r1.Start.J < r2.Start.J) return true;
            else if (r1.Start.J > r2.Start.J) return false;
            else
                if (r1.Start.K < r2.Start.K) return true;
                else if (r1.Start.K > r2.Start.K) return false;
                else
                    if (r1.End.I < r2.End.I) return true;
                    else if (r1.End.I > r2.End.I) return false;
                    else
                        if (r1.End.J < r2.End.J) return true;
                        else if (r1.End.J > r2.End.J) return false;
                        else
                            if (r1.End.K < r2.End.K) return true;
                            else if (r1.End.K > r2.End.K) return false;
                            else
                                return false;
    }

    std::ostream& Dump(std::ostream& o) const
    {
        o << BlockID() << ":" << MeshRange();
        return o;
    }

protected:
    friend class Roster;

    BlockPatch(int blockID, const IndexRange& meshRange, int uniqueID = -1)
    :   mBlockID(blockID), mMeshRange(meshRange), mCellRange(MeshRangeToCellRange(meshRange)), mUniqueID(uniqueID)
    {
        InitOrientation();
    }

    void InitOrientation();

    IndexRange MeshRangeToCellRange(const IndexRange& mr) const;

private:
    int mBlockID;
    IndexRange mMeshRange;
    IndexRange mCellRange;
    int mUniqueID;

    Direction mDir;
    IndexIJK mI1, mI2;
    IndexIJK mI3Ghost, mI3Interior;
    IndexIJK mDI3Ghost;
};

inline std::ostream& operator << (std::ostream& o, const BlockPatch& bp)
{
    return bp.Dump(o);
}

#endif // INCLUDED_BLOCK_PATCH_H__
