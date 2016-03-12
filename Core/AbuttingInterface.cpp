#include "Communicator.h"
#include "AbuttingInterface.h"
#include "Roster.h"
#include <algorithm>
#include <iterator>

AbuttingInterface::~AbuttingInterface()
{
    for (PatchMeshes::iterator i = mPatchMeshes.begin(); i != mPatchMeshes.end(); ++i)
    {
        Structured<double> pxyz = i->second;
        delete[] pxyz.Data;
    }
}

BlockPatches
AbuttingInterface::AllBlockPatches() const
{
    BlockPatches bps(mBlockPatches);
    copy(mDonorBlockPatches.begin(), mDonorBlockPatches.end(), std::back_inserter(bps));
    return bps;
}

std::set<int>
AbuttingInterface::BlockIDs() const
{
    BlockPatches bps = AllBlockPatches();

    std::set<int> blockIDs;
    for (BlockPatches::const_iterator i = bps.begin(); i != bps.end(); ++i)
    {
        std::cout << *i << std::endl;
        blockIDs.insert(i->BlockID());
    }

    return blockIDs;
}

void
AbuttingInterface::SetPatchMesh(int blockID, const Structured<double>& XYZ)
{
    std::ostream& o = Communicator::GetInstance()->Console();

    BlockPatches bps = AllBlockPatches();
    for (BlockPatches::const_iterator i = bps.begin(); i != bps.end(); ++i)
    {
        const BlockPatch& bp = *i;
        if (bp.BlockID() != blockID)
            continue;

        Structured<double>& pxyz = mPatchMeshes[bp.UniqueID()];
        pxyz.Allocate(3, bp.MeshRange());
        for (IndexIterator itor(bp.MeshRange()); !itor.IsEnd(); itor.Advance())
        {
            IndexIJK ijk = itor.Index();
            double* dst = pxyz(ijk);
            double* src = XYZ(ijk);
            dst[0] = src[0];
            dst[1] = src[1];
            dst[2] = src[2];
        }
        o << "AbbuttingInterface: read block " << blockID << " range " << bp.MeshRange() << std::endl;
    }
}

