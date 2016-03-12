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

