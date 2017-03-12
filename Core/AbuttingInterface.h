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
#ifndef INCLUDED_ABUTTING_INTERFACE_H__
#define INCLUDED_ABUTTING_INTERFACE_H__

#include "Structured.h"
#include "BlockPatch.h"
#include "Model.h"
#include "IterationContext.h"
#include <vector>
#include <set>
#include <map>
#include <iostream>

class InterfaceDataAdaptorBase;

class AbuttingInterface
{
public:
    typedef std::map<int, Structured<double> > PatchMeshes;

    virtual ~AbuttingInterface();

    virtual void InitializeMapping() = 0;
    virtual void MapMesh(const IterationContext& iteration) = 0;
    virtual void MapData(const Model& model, InterfaceDataAdaptorBase* adaptor) const = 0;
    virtual void DumpMappingResult(std::ostream& o) const = 0;

    const BlockPatches& SelfBlockPatches() const { return mBlockPatches; }
    const BlockPatches& DonorBlockPatches() const { return mDonorBlockPatches; }

    // returns all "self" and "donor" patches
    BlockPatches AllBlockPatches() const;

    std::set<int> BlockIDs() const;
    void SetPatchMesh(int blockID, const Structured<double>& XYZ); // FIXME: not elegant at all
    const PatchMeshes& GetPatchMeshes() const { return mPatchMeshes; }

protected:
    AbuttingInterface(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches)
    : mBlockPatches(blockPatches), mDonorBlockPatches(donorBlockPatches)
    {}

    PatchMeshes& GetPatchMeshes() { return mPatchMeshes; }
    void ConvertMappedDataToLocalFrame(const Model& model, Structured<double>& U, const BlockPatch& bp) const;
    void DumpInterfaceGhostCells(std::ostream& o, const Model& model) const;

private:
    BlockPatches mBlockPatches;
    BlockPatches mDonorBlockPatches;

    PatchMeshes mPatchMeshes;
};

typedef std::vector<AbuttingInterface*> Interfaces;

class InterfaceDataAdaptorBase
{
public:
    virtual ~InterfaceDataAdaptorBase() {}

    virtual Structured<double>& GetBlockData(int blockID) const = 0;
    virtual const Structured<double>& GetBlockPatchData(int blockPatchUniqueID) const = 0;

protected:

private:
};

#endif // INCLUDED_ABUTTING_INTERFACE_H__

