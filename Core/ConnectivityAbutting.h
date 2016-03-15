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
// $Id: ConnectivityAbutting.h 41 2010-07-06 10:48:07Z kato $
#ifndef INCLUDED_CONNECTIVITY_ABUTTING_H__
#define INCLUDED_CONNECTIVITY_ABUTTING_H__

#include "Connectivity.h"
#include "Connectivity1to1.h" // for IndexTransform
#include <string>
#include <vector>

#if 0
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
#endif

class ConnectivityAbutting : public Connectivity
{
public:
    ConnectivityAbutting(
        const IndexRange& meshRange, Direction dir, const Block& block,
        const IndexRange& donorMeshRange, int donorBlockID,
        int tag
        );
    virtual ~ConnectivityAbutting();

    virtual void Apply(const Block& block, Structured<double>& U) {}
    virtual void ApplyTurb(const Block& block, Structured<double>& UT, const Structured<double>& U) {}

#if 0
    void AddDonor(int donorBlockID, const char* name, const IndexRange& donorMeshRange)
    {
        mDonorPatches.push_back(DonorPatch(donorBlockID, name, donorMeshRange));
    }

    const DonorPatches& GetDonorPatches() const
    {
        return mDonorPatches;
    }
#endif

protected:

private:
    int mBlockID, mDonorBlockID;
    int mTag;
#if 0
    DonorPatches mDonorPatches;
#endif
};

#endif // INCLUDED_CONNECTIVITY_ABUTTING_H__

