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
// $Id: PlanarMapping.h 261 2013-01-05 12:01:20Z kato $
#ifndef INCLUDED_PLANAR_MAPPING_H__
#define INCLUDED_PLANAR_MAPPING_H__

#include "MeshMapper.h"
#include "BlockPatch.h"
#include "AbuttingInterface.h"
#include "Structured.h"
#include "Model.h"

class Vector3;

class PlanarMapping : public MeshMapper
{
public:
    PlanarMapping(const Structured<double>& patchXYZ, const BlockPatch& bp);
    virtual ~PlanarMapping();

    virtual void InitializeMapping();
    virtual void MapOn(const BlockPatch& bpd, const Structured<double>& XYZDonor, const Vector3& angleSelf, const Vector3& angleDonor, bool debug = false);
    virtual bool AllMapped() const;

    virtual void MapData(const Model& model, Structured<double>& U, const BlockPatch& bpd, const Structured<double>& UDonor) const;

    const Structured<int>& DonorCells() const { return mCells; }
    const Structured<double>& Centroids() const { return mC; }
    const Structured<double>& Distances() const { return mD; }

protected:
    Vector3 PatchNormal(const Structured<double>& XYZ) const;
    void PatchBasis(Vector3& e1, Vector3& e2, const Structured<double>& XYZ) const;
    //IndexRange CellFaceRange(const IndexRange& meshRange) const;
    //void CellFaceRange(IndexRange& cfr, IndexIJK& i1, IndexIJK& i2, IndexIJK& i3, const IndexRange& meshRange) const;
    void ComputeCentroids(const BlockPatch& bp, Structured<double>& centroids, const Structured<double>& XYZ) const;
    void InterpolateOn(double& c12, double& c13, const Vector3& p, const Structured<double>& XYZ,
            const IndexIJK& ijk, const IndexIJK& i1, const IndexIJK& i2) const;

private:
    const Structured<double>& mXYZ;
    BlockPatch mBP;
    Structured<int> mCells; // contains indices to the closest donor cells. for each face at [i, j, k], it has [DonorBlockID, iDonor, jDonor, kDonor]
    Structured<double> mC; // centroids
    Structured<double> mD; // distances
};

#endif // INCLUDED_PLANAR_MAPPING_H__

