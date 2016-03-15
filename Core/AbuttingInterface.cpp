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

