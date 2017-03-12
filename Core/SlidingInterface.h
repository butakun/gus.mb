/*
    gus.mb, an open source flow solver.
    Copyright (C) 2017 Hiromasa Kato <hiromasa at gmail.com>

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
#ifndef INCLUDED_SLIDING_INTERFACE_H__
#define INCLUDED_SLIDING_INTERFACE_H__

#include "AbuttingInterface.h"
#include "Structured.h"
#include "BlockPatch.h"
#include "Model.h"
#include "IterationContext.h"
#include "Roster.h"
#include <ANN/ANN.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <vector>
#include <iostream>

class InterfaceDataAdaptorBase;

class SlidingInterface : public AbuttingInterface
{
public:
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    //typedef CGAL::Simple_cartesian<double> Kernel;
    typedef Kernel::Point_2 Point_2;
    typedef Kernel::Segment_2 Segment_2;
    typedef CGAL::Polygon_2<Kernel> Polygon_2;

    class PatchCellFaceIndex
    {
    public:
        int BlockPatchID, I, J, K;
        PatchCellFaceIndex(int bpID, int i, int j, int k) : BlockPatchID(bpID), I(i), J(j), K(k) {}
    };

    class CellFaceIntersection
    {
    public:
        double Area;
        size_t CellFaceIndex;
        CellFaceIntersection(double area, size_t cellFaceIndex) : Area(area), CellFaceIndex(cellFaceIndex) {}
    };

    static SlidingInterface* New(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches);

    virtual ~SlidingInterface();

    virtual void InitializeMapping();
    virtual void MapMesh(const IterationContext& iteration);
    virtual void MapData(const Model& model, InterfaceDataAdaptorBase* adaptor) const;
    virtual void DumpMappingResult(std::ostream& o) const {}

protected:
    SlidingInterface(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches)
    : AbuttingInterface(blockPatches, donorBlockPatches), mDonorCentroids(NULL), mKDTree(NULL), mRadial(false)
    {}

    void DetectRadialInterface();

    Point_2 MapTo2DRadial(const Vector3& p) const; // deprecated
    void MapTo2DRadial(double& theta, double& eta, const Vector3& p) const;
    Polygon_2 MapTo2DRadial(int& quadrant, const Vector3* p) const; // quadrant = 1, 0, -1 (all theta positive, mixed, negative)
    void ComputeCentroid(ANNpoint centroid, const Polygon_2& poly) const;
    void DumpIntersection(size_t iSelf, ANNidxArray nnIdx, int k) const;

private:
    ANNpointArray mDonorCentroids;
    ANNkd_tree* mKDTree;

    std::vector<PatchCellFaceIndex> mSelfCellFaceIndices, mDonorCellFaceIndices;
    std::vector<Polygon_2> mSelfCellFaces, mDonorCellFaces;

    //std::vector<Kernel::FT> mSelfCellFaceAreas, mDonorCellFaceAreas;
    std::vector<double> mSelfCellFaceAreas, mDonorCellFaceAreas;
    std::vector<std::list<CellFaceIntersection> > mCellFaceIntersections;
    //std::vector<double> mSelfCellFaceSqLengths; // (squared) characteristic lengths used for ANN search

    bool mRadial;
    std::string mAxis;
    Vector3 mAxis1, mAxis2, mAxis3; // mAxis3 is the axis of rotation
    double mRadius;
};

inline
std::ostream& operator<< (std::ostream& o, const SlidingInterface::PatchCellFaceIndex& pcfi)
{
    o << Roster::GetInstance()->GetBlockPatch(pcfi.BlockPatchID).BlockID() << ":" << IndexIJK(pcfi.I, pcfi.J, pcfi.K);
    return o;
}

inline
std::ostream& operator<< (std::ostream& o, const SlidingInterface::Point_2& p)
{
    o << '(' << CGAL::to_double(p.x()) << ',' << CGAL::to_double(p.y()) << ')';
    return o;
}

inline
std::ostream& operator<< (std::ostream& o, const SlidingInterface::Polygon_2& poly)
{
    o << poly[0] << ',' << poly[1] << ',' << poly[2] << ',' << poly[3];
    return o;
}


#endif // INCLUDED_SLIDING_INTERFACE_H__

