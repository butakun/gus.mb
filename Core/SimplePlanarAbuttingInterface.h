#ifndef INCLUDED_SIMPLE_PLANAR_ABUTTING_INTERFACE_H__
#define INCLUDED_SIMPLE_PLANAR_ABUTTING_INTERFACE_H__

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

class SimplePlanarAbuttingInterface
{
public:
    typedef std::vector<BlockPatch> BlockPatches;
    typedef std::map<int, Structured<double> > PatchMeshes;
    typedef std::map<int, PlanarMapping*> Mappings;

    static SimplePlanarAbuttingInterface* New(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches);

    virtual ~SimplePlanarAbuttingInterface();

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
    SimplePlanarAbuttingInterface(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches)
    : mBlockPatches(blockPatches), mDonorBlockPatches(donorBlockPatches)
    {}

    void ConvertMappedDataToLocalFrame(const Model& model, Structured<double>& U, const BlockPatch& bp) const;

private:
    BlockPatches mBlockPatches;
    BlockPatches mDonorBlockPatches;

    PatchMeshes mPatchMeshes;

    Mappings mMappings;
};

typedef std::vector<SimplePlanarAbuttingInterface*> Interfaces;

class InterfaceDataAdaptorBase
{
public:
    virtual ~InterfaceDataAdaptorBase() {}

    virtual Structured<double>& GetBlockData(int blockID) const = 0;
    virtual const Structured<double>& GetBlockPatchData(int blockPatchUniqueID) const = 0;

protected:

private:
};

#endif // INCLUDED_SIMPLE_PLANAR_ABUTTING_INTERFACE_H__

