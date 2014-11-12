// $Id: ConnectivityAbutting.h 41 2010-07-06 10:48:07Z kato $
#ifndef INCLUDED_CONNECTIVITY_ABUTTING_H__
#define INCLUDED_CONNECTIVITY_ABUTTING_H__

#include "Connectivity.h"
#include <string>
#include <vector>

class DonorPatch
{
public:
    int DonorBlockID;
    std::string Name;
    IndexRange MeshRange;

    DonorPatch(int blockID, const char* name, const IndexRange& meshRange)
    : DonorBlockID(blockID), Name(name), MeshRange(meshRange) {}
};

typedef std::vector<DonorPatch> DonorPatches;

class ConnectivityAbutting : public Connectivity
{
public:
    ConnectivityAbutting(const IndexRange& meshRange, Direction dir)
    : Connectivity(meshRange, dir)
    {}

    virtual ~ConnectivityAbutting()
    {}

    void AddDonor(int donorBlockID, const char* name, const IndexRange& donorMeshRange);
    {
        mDonorPatches.push_back(DonorPatch(donorBlockID, name, donorMeshRange));
    }

    const DonorPatches& GetDonorPatches() const
    {
        return mDonorPatches;
    }

protected:

private:
    DonorPatches mDonorPatches;
};

#endif // INCLUDED_CONNECTIVITY_ABUTTING_H__

