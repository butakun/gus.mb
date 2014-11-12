// $Id: AbuttingInterface.h 261 2013-01-05 12:01:20Z kato $
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

class PlanarMapping;
class InterfaceDataAdaptorBase;

class AbuttingInterface
{
public:
    typedef std::vector<BlockPatch> BlockPatches;
    typedef std::map<int, Structured<double> > PatchMeshes;
    typedef std::map<int, PlanarMapping*> Mappings;

    static AbuttingInterface* New(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches);

    virtual ~AbuttingInterface();

    const BlockPatches& SelfBlockPatches() const { return mBlockPatches; }
    const BlockPatches& DonorBlockPatches() const { return mDonorBlockPatches; }

    // returns all "self" and "donor" patches
    BlockPatches AllBlockPatches() const;

    std::set<int> BlockIDs() const;
    void SetPatchMesh(int blockID, const Structured<double>& XYZ); // FIXME: not elegant at all
    const PatchMeshes& GetPatchMeshes() const { return mPatchMeshes; }

    void InitializeMapping();
    void MapMesh(const IterationContext& iteration);
    void MapData(const Model& model, InterfaceDataAdaptorBase* adaptor) const;
    const PlanarMapping& GetMapper(int blockPatchUniqueID) const;

    void DumpMappingResult(std::ostream& o) const;

protected:
    AbuttingInterface(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches)
    : mBlockPatches(blockPatches), mDonorBlockPatches(donorBlockPatches)
    {}

    void ConvertMappedDataToLocalFrame(const Model& model, Structured<double>& U, const BlockPatch& bp) const;

private:
    BlockPatches mBlockPatches;
    BlockPatches mDonorBlockPatches;

    PatchMeshes mPatchMeshes;

    Mappings mMappings;
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

