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
#ifndef INCLUDED_BITMAP_MESH_INTERPOLATOR_H__
#define INCLUDED_BITMAP_MESH_INTERPOLATOR_H__

#include "BlockPatch.h"
#include "Vector3.h"
#include <vector>
#include <map>

class CellLocator
{
public:
    typedef std::vector<BlockPatch> BlockPatches;
    typedef std::map<int, Structured<double> > PatchMeshMap;

    struct PatchCell {
        int id;
        IndexIJK index;
        PatchCell(int id_, const IndexIJK& index_) : id(id_), index(index_) {}
        PatchCell() : id(-1) {}
    };

    CellLocator() {}
    virtual ~CellLocator() {}

    void AddPatch(const BlockPatch& bp, const Structured<double>& XYZ);
    bool LocatePoint(PatchCell& cell, const Vector3& p) const;

    bool LocatePointInPatch(IndexIJK& ijk, const Vector3& p, const BlockPatch& bp, const Structured<double>& XYZ) const;

    // p4 --- p3
    // |   p   |
    // p1 --- p2
    bool IsPointInside(const Vector3& p, const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4) const;

    //    p3
    //   /p|
    // p1--p2
    bool IsPointInsideTriangle(const Vector3& p, const Vector3& p1, const Vector3& p2, const Vector3& p3) const;

protected:

private:
    BlockPatches mBlockPatches;
    PatchMeshMap mPatchMeshMap;
};

class MeshInterpolator
{
public:
    virtual ~MeshInterpolator() {}
};

class BitmapMeshInterpolator : public MeshInterpolator
{
public:
    typedef std::vector<BlockPatch> BlockPatches;
    typedef std::map<int, Structured<double> > PatchMeshes;

    class PatchIndex {
    public:
        int PatchID[2];
        IndexIJK PatchIJK[2];

        PatchIndex() {}
        PatchIndex(const PatchIndex& pi)
        {
            PatchID[0] = pi.PatchID[0]; PatchID[1] = pi.PatchID[1];
            PatchIJK[0] = pi.PatchIJK[0]; PatchIJK[1] = pi.PatchIJK[1];
        }
    };

    class DonorCellFace {
    public:
        int PatchID;
        IndexIJK Index;
        double Weight;
        DonorCellFace(int id, const IndexIJK& ijk, double w) : PatchID(id), Index(ijk), Weight(w) {}
    };
    typedef std::vector<DonorCellFace> DonorCellFaces;

    BitmapMeshInterpolator();
    virtual ~BitmapMeshInterpolator();

    void AddPatch(int which, const BlockPatch& bp, const Structured<double>& XYZ); // which = 0 (mPatches1) or 1 (mPatches2)

//protected:
    void GenerateBitmapPixels();
    void ComputeCellFaceWeights();
    void ComputeMappingWeights();
    const std::vector<Vector3>& GetPixels() const { return mPixels; }
    const std::vector<PatchIndex>& GetPixelMapping() const { return mPixelMapping; }
    std::map<int, Structured<double> >& GetCellFaceWeights() { return mCellFaceWeights; }

private:
    BlockPatches mPatches1, mPatches2;
    PatchMeshes mPatchMeshes;

    std::vector<Vector3> mPixels;           // each point represents a pixel
    std::vector<double> mPixelWeights;      // pixels are not uniformly distributed, so each point has weight assigned to it that is proportional to the pixel area
    std::vector<PatchIndex> mPixelMapping;  // for each pixel, stores the indices to the cells on Patch1 and Patch2 that contain the pixel

    std::map<int, Structured<double> > mCellFaceWeights;    // sum of the pixels mapped to each cell
    std::map<int, Structured<size_t> > mCellFaceMapping;
    std::vector<DonorCellFaces> mDonorList;
};

#endif // INCLUDED_BITMAP_MESH_INTERPOLATOR_H__

