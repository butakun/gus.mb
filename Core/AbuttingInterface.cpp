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

void
AbuttingInterface::ConvertMappedDataToLocalFrame(const Model& model, Structured<double>& U, const BlockPatch& bp) const
{
    //std::ostream& LOG = Communicator::GetInstance()->Console();

    const Block& block = dynamic_cast<const Block&>(*Roster::GetInstance()->GetBlock(bp.BlockID()));

    IndexRange cr0;
    IndexIJK i1, i2, di3;
    bp.GhostCellRange(cr0, i1, i2, di3);

    int ndof = U.DOF();
    double * tmp = new double[ndof];

    for (IndexIterator itor(cr0); !itor.IsEnd(); itor.Advance())
    {
        IndexIJK ijk = itor.Index();
        double* u = U(ijk);
        for (int l = 0; l < ndof; ++l)
        {
            tmp[l] = u[l];
        }
        model.FromGlobalToLocal(u, tmp, block, ijk);
    }
    for (int irind = 1; irind < block.GhostLayers(); ++irind)
    {
        for (IndexIterator itor(cr0); !itor.IsEnd(); itor.Advance())
        {
            IndexIJK ijk = itor.Index();
            double* u0 = U(ijk);
            double* u = U(ijk + di3 * irind);
            for (int l = 0; l < ndof; ++l)
            {
                u[l] = u0[l];
            }
        }
    }

    delete[] tmp;
}

void
AbuttingInterface::DumpInterfaceGhostCells(std::ostream& o, const Model& model) const
{
    for (const auto& bp : SelfBlockPatches())
    {
        o << "DumpInterfaceGhostCells: " << bp << std::endl;
        const Block& block = dynamic_cast<const Block&>(*Roster::GetInstance()->GetBlock(bp.BlockID()));
        const Structured<double>& U = block.U(); // FIXME
        size_t dof = U.DOF();

        IndexRange cr;
        IndexIJK i1, i2, di3;
        bp.GhostCellRange(cr, i1, i2, di3);

        for (IndexIterator itor(cr); !itor.IsEnd(); itor.Advance())
        {
            IndexIJK ijk = itor.Index();
            double* u = U(ijk);
            o << ijk << ':';
            for (int l = 0; l < dof; ++l)
            {
                o << u[l] << ' ';
            }
            o << std::endl;
        }
    }
}

