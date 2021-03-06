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
#include "SimplePlanarAbuttingInterface.h"
#include "PlanarMapping.h"
#include "Roster.h"
#include "RigidBodyMotion.h"
#include <algorithm>
#include <iterator>

SimplePlanarAbuttingInterface*
SimplePlanarAbuttingInterface::New(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches)
{
    SimplePlanarAbuttingInterface* interface = new SimplePlanarAbuttingInterface(blockPatches, donorBlockPatches);
    Roster::GetInstance()->AddInterface(interface);

    return interface;
}

SimplePlanarAbuttingInterface::~SimplePlanarAbuttingInterface()
{
#if 0
    for (PatchMeshes::iterator i = mPatchMeshes.begin(); i != mPatchMeshes.end(); ++i)
    {
        Structured<double> pxyz = i->second;
        delete[] pxyz.Data;
    }
#endif

    for (Mappings::iterator i = mMappings.begin(); i != mMappings.end(); ++i)
    {
        PlanarMapping* mapping = i->second;
        delete mapping;
    }
}

#if 0
SimplePlanarAbuttingInterface::BlockPatches
SimplePlanarAbuttingInterface::AllBlockPatches() const
{
    BlockPatches bps(mBlockPatches);
    copy(mDonorBlockPatches.begin(), mDonorBlockPatches.end(), std::back_inserter(bps));
    return bps;
}

std::set<int>
SimplePlanarAbuttingInterface::BlockIDs() const
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
SimplePlanarAbuttingInterface::SetPatchMesh(int blockID, const Structured<double>& XYZ)
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
#endif

void
SimplePlanarAbuttingInterface::InitializeMapping()
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    for (BlockPatches::const_iterator i = SelfBlockPatches().begin();
        i != SelfBlockPatches().end(); ++i)
    {
        const BlockPatch& bp = *i;
        PatchMeshes::const_iterator ipm = GetPatchMeshes().find(bp.UniqueID());
        assert(ipm != GetPatchMeshes().end());
        const Structured<double>& XYZ = ipm->second;

        // Make sure there is no mapper created already for this patch.
        Mappings::const_iterator im = mMappings.find(bp.UniqueID());
        assert(im == mMappings.end());

        // Store the mapper, indexed with patch's unique id.
        PlanarMapping* mapping = new PlanarMapping(XYZ, bp);
        LOG << "AbuttinInterface::InitializeMapping: created a mapper for " << bp << " UniqueID = " << bp.UniqueID() << std::endl;
        mMappings[bp.UniqueID()] = mapping;
    }
}

void
SimplePlanarAbuttingInterface::MapMesh(const IterationContext& iteration)
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    LOG << "SimplePlanarAbuttingInterface::MapMesh" << std::endl;

    // Initialize the mapping information associated with the patches of this interface.
    for (Mappings::iterator i = mMappings.begin();
        i != mMappings.end(); ++i)
    {
        PlanarMapping* mapping = i->second;
        mapping->InitializeMapping();
    }

    // Visit all the patches one by one, and map against the donor patches.
    for (BlockPatches::const_iterator i = SelfBlockPatches().begin();
        i != SelfBlockPatches().end(); ++i)
    {
        const BlockPatch& bp = *i;
        Mappings::iterator im = mMappings.find(bp.UniqueID());
        assert(im != mMappings.end());
        PlanarMapping* mapping = im->second;

        // Block rotation angle at this instant
        const VirtualBlock* block = Roster::GetInstance()->GetBlock(bp.BlockID());
        Vector3 angle0(0.0, 0.0, 0.0);
        if (!block->IsStationary())
        {
            RotationalMotion *rm = dynamic_cast<RotationalMotion*>(block->GetRigidBodyMotion());
            assert(rm); // FIXME: what about translational motion?
            angle0 = rm->AngularVelocity() * iteration.Time();
            LOG << "SimplePlanarAbuttingInterface::MapMesh: I am a rotating domain, currently at angle0 = " << angle0 * 180.0 / M_PI << " deg., patch = " << bp << std::endl;
        }

        // We keep track of how much each block has been rotated.
        // If rotated more than 360deg, we won't bother with that block.
        std::map<int, int> blockAngles;
        for (BlockPatches::const_iterator j = DonorBlockPatches().begin();
            j != DonorBlockPatches().end(); ++j)
        {
            blockAngles[j->BlockID()] = 0;
        }

        bool allMapped = false;
        //while (!allMapped)
        for ( ; ; )
        {
            size_t countBlocksDone = 0;
            for (BlockPatches::const_iterator j = DonorBlockPatches().begin();
                j != DonorBlockPatches().end(); ++j)
            {
                const BlockPatch& bpd = *j;
                const VirtualBlock* donorBlock = Roster::GetInstance()->GetBlock(bpd.BlockID());
                const Structured<double>& XYZD = GetPatchMeshes()[bpd.UniqueID()];

                // Donor block rotation
                Vector3 angle0d(0.0, 0.0, 0.0);
                if (!donorBlock->IsStationary())
                {
                    RotationalMotion *rmd = dynamic_cast<RotationalMotion*>(donorBlock->GetRigidBodyMotion());
                    assert(rmd); // FIXME: what about translational motion?
                    angle0d = rmd->AngularVelocity() * iteration.Time();
                    LOG << "SimplePlanarAbuttingInterface::MapMesh: donor is a rotating domain, currently at angle0d = " << angle0d * 180.0 / M_PI << " deg., patch = " << bpd << std::endl;
                }

                Vector3 angled(0.0, 0.0, 0.0);
                if (donorBlock->IsPeriodic())
                {
                    angled = donorBlock->Periodicity() * blockAngles[bpd.BlockID()];
                    LOG << "Rotating block " << bpd.BlockID() << " by " << angled * 180.0 / M_PI << std::endl;
                    if (angled.Mag() >= 2.0 * M_PI)
                    {
                        countBlocksDone++;
                        LOG << "Skipping block " << bpd.BlockID() << " (been rotated " << angled.Mag() * 180.0 / M_PI << " deg." << std::endl;
                        continue;
                    }
                }
                else
                {
                    countBlocksDone++;
                }

                mapping->MapOn(bpd, XYZD, angle0, angle0d + angled, false);
                allMapped = mapping->AllMapped();
                blockAngles[bpd.BlockID()] += 1;
            }
            if (countBlocksDone == DonorBlockPatches().size())
            {
                if (!allMapped)
                {
                    LOG << "SimplePlanarAbuttingInterface: *** some interface patches are not completely matched. ***" << std::endl;
                }
                break;
            }
        }
    }
}

void
SimplePlanarAbuttingInterface::MapData(const Model& model, InterfaceDataAdaptorBase* adaptor) const
{
    for (BlockPatches::const_iterator i = SelfBlockPatches().begin();
        i != SelfBlockPatches().end(); ++i)
    {
        const BlockPatch& bp = *i;
        Structured<double>& U = adaptor->GetBlockData(bp.BlockID());
        const PlanarMapping* mapping = mMappings.find(bp.UniqueID())->second;;
        assert(mapping != NULL);

        for (BlockPatches::const_iterator j = DonorBlockPatches().begin();
            j != DonorBlockPatches().end(); ++j)
        {
            const BlockPatch& dbp = *j;
            const Structured<double>& Ud = adaptor->GetBlockPatchData(dbp.UniqueID());
            mapping->MapData(model, U, dbp, Ud);
        }

        ConvertMappedDataToLocalFrame(model, U, bp);
    }
}

const PlanarMapping&
SimplePlanarAbuttingInterface::GetMapper(int blockPatchUniqueID) const
{
    Mappings::const_iterator i = mMappings.find(blockPatchUniqueID);
    assert(i != mMappings.end());
    return *(i->second);
}

void
SimplePlanarAbuttingInterface::DumpMappingResult(std::ostream& o) const
{
    o << "<VTKFile type=\"PolyData\">" << std::endl;
    o << "  <PolyData>" << std::endl;
    for (BlockPatches::const_iterator ibp = SelfBlockPatches().begin(); ibp != SelfBlockPatches().end(); ++ibp)
    {
        const BlockPatch& bp = *ibp;
        const PlanarMapping* mapping = mMappings.find(bp.UniqueID())->second;;
        const IndexRange& cfr = bp.CellFaceRange();

        o << "    <Piece NumberOfPoints=\"" << cfr.Count() << "\">" << std::endl;
        o << "      <!-- Block " << bp.BlockID() << " MeshRange = " << bp.MeshRange() << " -->" << std::endl;
        o << "      <Points>" << std::endl;
        o << "        <DataArray type = \"Float32\" format=\"ascii\" NumberOfComponents=\"3\">" << std::endl;
        const Structured<int>& donorCells = mapping->DonorCells();
        const Structured<double>& centroids = mapping->Centroids();
        for (IndexIterator ip(cfr); !ip.IsEnd(); ip.Advance())
        {
            const IndexIJK& ijk = ip.Index();
            double* c = centroids(ijk);
            o << c[0] << ' ' << c[1] << ' ' << c[2] << std::endl;
        }
        o << "        </DataArray>" << std::endl;
        o << "      </Points>" << std::endl;
        o << "      <PointData>" << std::endl;
        o << "        <DataArray type=\"UInt32\" format=\"ascii\" Name=\"DonorBlockID\">" << std::endl;
        for (IndexIterator ip(cfr); !ip.IsEnd(); ip.Advance())
        {
            const IndexIJK& ijk = ip.Index();
            int* d = donorCells(ijk);
            o << d[0] << std::endl;
        }
        o << "        </DataArray>" << std::endl;
        o << "      </PointData>" << std::endl;
        o << "    </Piece>" << std::endl;
    }
    o << "  </PolyData>" << std::endl;
    o << "</VTKFile>" << std::endl;
}

