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
#ifndef INCLUDED_SIMPLE_PLANAR_ABUTTING_INTERFACE_H__
#define INCLUDED_SIMPLE_PLANAR_ABUTTING_INTERFACE_H__

#include "AbuttingInterface.h"
#include "Structured.h"
#include "BlockPatch.h"
#include "Model.h"
#include "IterationContext.h"
#include <vector>
#include <set>
#include <map>
#include <iostream>

class PlanarMapping;
class InterfaceDataAdaptorBase;

class SimplePlanarAbuttingInterface : public AbuttingInterface
{
public:
    typedef std::map<int, PlanarMapping*> Mappings;

    static SimplePlanarAbuttingInterface* New(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches);

    virtual ~SimplePlanarAbuttingInterface();

    virtual void InitializeMapping();
    virtual void MapMesh(const IterationContext& iteration);
    virtual void MapData(const Model& model, InterfaceDataAdaptorBase* adaptor) const;

#if 0
    const BlockPatches& SelfBlockPatches() const { return mBlockPatches; }
    const BlockPatches& DonorBlockPatches() const { return mDonorBlockPatches; }

    // returns all "self" and "donor" patches
    BlockPatches AllBlockPatches() const;

    std::set<int> BlockIDs() const;
    void SetPatchMesh(int blockID, const Structured<double>& XYZ); // FIXME: not elegant at all
    const PatchMeshes& GetPatchMeshes() const { return mPatchMeshes; }
#endif

    const PlanarMapping& GetMapper(int blockPatchUniqueID) const;

    void DumpMappingResult(std::ostream& o) const;

protected:
    SimplePlanarAbuttingInterface(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches)
    : AbuttingInterface(blockPatches, donorBlockPatches)
    {}

    void ConvertMappedDataToLocalFrame(const Model& model, Structured<double>& U, const BlockPatch& bp) const;

private:
    Mappings mMappings;
};

#endif // INCLUDED_SIMPLE_PLANAR_ABUTTING_INTERFACE_H__

