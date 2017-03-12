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

#include "SlidingInterface.h"
#include "Roster.h"
#include "Communicator.h"
#include <boost/range/combine.hpp>

#define SLIDING_INTERFACE_DEBUG 0

SlidingInterface*
SlidingInterface::New(const BlockPatches& blockPatches, const BlockPatches& donorBlockPatches)
{
    SlidingInterface* interface = new SlidingInterface(blockPatches, donorBlockPatches);
    Roster::GetInstance()->AddInterface(interface);
    return interface;
}

SlidingInterface::~SlidingInterface()
{
    delete mKDTree;
    annDeallocPts(mDonorCentroids);
}

void
SlidingInterface::DetectRadialInterface()
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    double varRx, varRy, varRz; // variance in Rx, Ry, Rz.
    double avgRx, avgRy, avgRz;

    // Gather all the patches, both local and remote.
    size_t numVertices = 0;
    std::vector<BlockPatch> allbps;
    for (BlockPatches::const_iterator i = SelfBlockPatches().begin();
        i != SelfBlockPatches().end(); ++i)
    {
        numVertices += i->MeshRange().Count();
        allbps.push_back(*i);
    }
    for (BlockPatches::const_iterator i = DonorBlockPatches().begin();
        i != DonorBlockPatches().end(); ++i)
    {
        numVertices += i->MeshRange().Count();
        allbps.push_back(*i);
    }

    double* Rs = new double[numVertices * 3];
    double* rs = Rs;
    avgRx = avgRy = avgRz = 0.0;
    for (std::vector<BlockPatch>::const_iterator i = allbps.begin();
        i != allbps.end(); ++i)
    {
        const BlockPatch& bp = *i;
        const Structured<double>& XYZ = GetPatchMeshes()[bp.UniqueID()];
        IndexRange mr = bp.MeshRange();

        for (IndexIterator mritor(mr); !mritor.IsEnd(); mritor.Advance())
        {
            Vector3 p(XYZ(mritor.Index()));
            double Rx, Ry, Rz;
            Rx = std::sqrt(p.Y() * p.Y() + p.Z() * p.Z());
            Ry = std::sqrt(p.X() * p.X() + p.Z() * p.Z());
            Rz = std::sqrt(p.X() * p.X() + p.Y() * p.Y());
            *rs++ = Rx;
            *rs++ = Ry;
            *rs++ = Rz;
            avgRx += Rx;
            avgRy += Ry;
            avgRz += Rz;
        }
    }
    avgRx /= double(numVertices);
    avgRy /= double(numVertices);
    avgRz /= double(numVertices);

    varRx = varRy = varRz = 0.0;
    rs = Rs;
    for (size_t i = 0; i < numVertices; ++i)
    {
        double drx, dry, drz;
        drx = *rs++ - avgRx;
        dry = *rs++ - avgRy;
        drz = *rs++ - avgRz;
        varRx += drx * drx;
        varRy += dry * dry;
        varRz += drz * drz;
    }
    varRx /= double(numVertices);
    varRy /= double(numVertices);
    varRz /= double(numVertices);

    LOG << "Detecting radial sliding mesh:"
        << " Var(Rx) = " << varRx
        << " Var(Ry) = " << varRy
        << " Var(Rz) = " << varRz
        << " std(Rx)/avg(Rx) = " << std::sqrt(varRx) / avgRx
        << " std(Ry)/avg(Ry) = " << std::sqrt(varRy) / avgRy
        << " std(Rz)/avg(Rz) = " << std::sqrt(varRz) / avgRz
        << std::endl;

    bool isAxisX, isAxisY, isAxisZ;
    isAxisX = std::sqrt(varRx) / avgRx < 1e-10;
    isAxisY = std::sqrt(varRy) / avgRy < 1e-10;
    isAxisZ = std::sqrt(varRz) / avgRz < 1e-10;
    if (isAxisX && !isAxisY && !isAxisZ)
    {
        mRadial = true;
        mAxis1 = Vector3(0, 1, 0);
        mAxis2 = Vector3(0, 0, 1);
        mAxis3 = Vector3(1, 0, 0);
        mRadius = avgRx;
    }
    if (!isAxisX && isAxisY && !isAxisZ)
    {
        mRadial = true;
        mAxis1 = Vector3(1, 0, 0);
        mAxis2 = Vector3(0, 0, 1);
        mAxis3 = Vector3(0, 1, 0);
        mRadius = avgRy;
    }
    if (!isAxisX && !isAxisY && isAxisZ)
    {   
        mRadial = true;
        mAxis1 = Vector3(1, 0, 0);
        mAxis2 = Vector3(0, 1, 0);
        mAxis3 = Vector3(0, 0, 1);
        mRadius = avgRz;
    }

    if (mRadial)
    {
        LOG << "Detected Axis1 = " << mAxis1 << ", Axis2 = " << mAxis2 << ", mAxis3 = " << mAxis3 << ", Radius = " << mRadius << std::endl;
    }
    else
    {
        LOG << "Not a radial interface" << std::endl;
    }

    delete[] Rs;
}

SlidingInterface::Point_2
SlidingInterface::MapTo2DRadial(const Vector3& p) const
{
    // xi = r * theta, eta = axial coordinate
    double eta = dot_product(p, mAxis3);
    double x1 = dot_product(p, mAxis1);
    double x2 = dot_product(p, mAxis2);
    double theta = std::atan2(x2, x1);
    double xi = mRadius * theta;
    return Point_2(xi, eta);
}

void
SlidingInterface::MapTo2DRadial(double& theta, double& eta, const Vector3& p) const
{
    eta = dot_product(p, mAxis3);
    double x1 = dot_product(p, mAxis1);
    double x2 = dot_product(p, mAxis2);
    theta = std::atan2(x2, x1);
}

SlidingInterface::Polygon_2
SlidingInterface::MapTo2DRadial(int& quadrant, const Vector3* p) const
{
    double theta[4], eta[4];
    Point_2 pts[4];
    for (int i = 0; i < 4; ++i)
    {
        MapTo2DRadial(theta[i], eta[i], p[i]);
    }
    if (theta[0] >= 0.0 && theta[1] >= 0.0 && theta[2] >= 0.0 && theta[3] >= 0.0)
    {
        quadrant = 1;
    }
    else if (theta[0] < 0.0 && theta[1] < 0.0 && theta[2] < 0.0 && theta[3] < 0.0)
    {
        quadrant = -1;
        for (int i = 0; i < 4; ++i)
            theta[i] = 2.0 * M_PI + theta[i];
    }
    else
    {
        quadrant = 0;
        for (int i = 0; i < 4; ++i)
        {
            if (theta[i] < -M_PI_2)
                theta[i] = 2.0 * M_PI + theta[i];
        }
    }

    for (int i = 0; i < 4; ++i)
    {
        pts[i] = Point_2(mRadius * theta[i], eta[i]);
    }
    return Polygon_2(pts, pts + 4);
}

void
SlidingInterface::ComputeCentroid(ANNpoint centroid, const Polygon_2& poly) const
{
    double p1[] = { CGAL::to_double(poly[0].x()), CGAL::to_double(poly[0].y()) };
    double p2[] = { CGAL::to_double(poly[1].x()), CGAL::to_double(poly[1].y()) };
    double p3[] = { CGAL::to_double(poly[2].x()), CGAL::to_double(poly[2].y()) };
    double p4[] = { CGAL::to_double(poly[3].x()), CGAL::to_double(poly[3].y()) };
    centroid[0] = 0.25 * (p1[0] + p2[0] + p3[0] + p4[0]);
    centroid[1] = 0.25 * (p1[1] + p2[1] + p3[1] + p4[1]);
}

void
SlidingInterface::InitializeMapping()
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    DetectRadialInterface();
    if (!mRadial) assert(false);

    std::vector<const BlockPatches*> bps_ = { &SelfBlockPatches(), &DonorBlockPatches() };
    std::vector<std::vector<PatchCellFaceIndex>*> cfis_ = { &mSelfCellFaceIndices, &mDonorCellFaceIndices };
    std::vector<std::vector<Polygon_2>*> cfs_ = { &mSelfCellFaces, &mDonorCellFaces };
    //std::vector<std::vector<Kernel::FT>*> cfas_ = { &mSelfCellFaceAreas, &mDonorCellFaceAreas };
    std::vector<std::vector<double>*> cfas_ = { &mSelfCellFaceAreas, &mDonorCellFaceAreas };

    for (auto tup : boost::combine(bps_, cfis_, cfs_, cfas_))
    {
        const BlockPatches* bps;
        std::vector<PatchCellFaceIndex>* cfis;
        std::vector<Polygon_2>* cfs;
        //std::vector<Kernel::FT>* cfas;
        std::vector<double>* cfas;
        boost::tie(bps, cfis, cfs, cfas) = tup;
        for (auto bp : *bps)
        {
            PatchMeshes::const_iterator ipm = GetPatchMeshes().find(bp.UniqueID());
            const Structured<double>& XYZ = ipm->second;

            IndexRange cfr = bp.CellFaceRange();
            IndexIJK i1 = bp.I1();
            IndexIJK i2 = bp.I2();
            for (IndexIterator cfitor(cfr); !cfitor.IsEnd(); cfitor.Advance())
            {
                IndexIJK ijk = cfitor.Index();
                Vector3 p[4];
                p[0] = Vector3(XYZ(ijk - i1 - i2));
                p[1] = Vector3(XYZ(ijk - i2));
                p[2] = Vector3(XYZ(ijk));
                p[3] = Vector3(XYZ(ijk - i1));

                int quadrant = 0;
                Polygon_2 poly = MapTo2DRadial(quadrant, p);
                cfs->push_back(poly);
                cfis->push_back(PatchCellFaceIndex(bp.UniqueID(), ijk.I, ijk.J, ijk.K));
                cfas->push_back(std::abs(CGAL::to_double(poly.area())));
            }
        }
    }
    LOG << "PatchCellFaceIndices: " << mSelfCellFaceIndices.size() << ", " << mDonorCellFaceIndices.size() << std::endl;
    LOG << "CellFaces: " << mSelfCellFaces.size() << ", " << mDonorCellFaces.size() << std::endl;
    LOG << "CellFaceAreas: " << mSelfCellFaceAreas.size() << ", " << mDonorCellFaceAreas.size() << std::endl;

    // add copies of donor cells beyond theta = [-pi:pi] FIXME
    size_t donorCellFacesCount = mDonorCellFaces.size();
    for (size_t i = 0; i < donorCellFacesCount; ++i)
    {
        const Polygon_2& polyOrig = mDonorCellFaces[i];
        double theta[4];
        for (size_t l = 0; l < 4; ++l)
            theta[l] = CGAL::to_double(polyOrig[l].x()) / mRadius;

        Point_2 pts[4];
        if (theta[0] > 1.5 * M_PI || theta[1] > 1.5 * M_PI || theta[2] > 1.5 * M_PI || theta[3] > 1.5 * M_PI)
        {
            for (size_t l = 0; l < 4; ++l)
                pts[l] = Point_2(mRadius * (theta[l] - 2.0 * M_PI), polyOrig[l].y());
        }
        else if (theta[0] < M_PI_2 || theta[1] < M_PI_2 || theta[2] < M_PI_2 || theta[3] < M_PI_2)
        {
            for (size_t l = 0; l < 4; ++l)
                pts[l] = Point_2(mRadius * (theta[l] + 2.0 * M_PI), polyOrig[l].y());
        }
        else
        {
            continue;
        }

        Polygon_2 poly(pts, pts + 4);
        //LOG << "Adding a copy of polygon " << i << " cyclically " << poly << std::endl;
        mDonorCellFaces.push_back(poly);
        mDonorCellFaceIndices.push_back(mDonorCellFaceIndices[i]);
        mDonorCellFaceAreas.push_back(mDonorCellFaceAreas[i]);
    }

    // build ANN database
#if 0
    mSelfCellFaceSqLengths.resize(mSelfCellFaces.size());
    for (size_t i = 0; i < mSelfCellFaces.size(); ++i)
    {
        const Polygon_2& poly = mSelfCellFaces[i];
        double sqLength = 0.0;
        for (Polygon_2::Edge_const_iterator j = poly.edges_begin();
            j != poly.edges_end(); ++j)
        {
            sqLength = std::max(sqLength, CGAL::to_double(j->squared_length()));
        }
        LOG << "Self cell " << i << " sqLength = " << sqLength << std::endl;
        mSelfCellFaceSqLengths[i] = sqLength;
    }
#endif

    mDonorCentroids = annAllocPts(mDonorCellFaces.size(), 2);
    for (size_t i = 0; i < mDonorCellFaces.size(); ++i)
    {
        const Polygon_2& poly = mDonorCellFaces[i];
        ComputeCentroid(mDonorCentroids[i], poly);
    }
    mKDTree = new ANNkd_tree(mDonorCentroids, mDonorCellFaces.size(), 2);

    // find intersections
    typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;

    ANNpoint selfCentroid = annAllocPt(2);
    double eps = 0.0;

    int kMax = 1000;
    ANNidxArray nnIdx = new ANNidx[kMax];
    ANNdistArray dists = new ANNdist[kMax];

    mCellFaceIntersections.resize(mSelfCellFaces.size());
    for (size_t i = 0; i < mSelfCellFaces.size(); ++i)
    {
        const Polygon_2& selfPoly = mSelfCellFaces[i];
        ComputeCentroid(selfCentroid, selfPoly);
        double selfArea = mSelfCellFaceAreas[i];
        //LOG << "Cell " << i << "(" << selfArea << ")" << selfPoly << std::endl;

        double areaIntersected = 0.0, overlap = 0.0, overlapMin = 0.95;

        int k = 0, dk = 10, kStart = 0;
        while (overlap < overlapMin)
        {
            k += dk;
            if (k > kMax)
            {
                DumpIntersection(i, nnIdx, kMax);
                exit(EXIT_FAILURE);
            }
            //LOG << " searching for " << k << " nearest neighbors" << std::endl;
            mKDTree->annkSearch(selfCentroid, k, nnIdx, dists, eps);

            for (size_t j = kStart; j < k; ++j)
            {
                size_t iDonor = nnIdx[j];
                const Polygon_2& donorPoly = mDonorCellFaces[iDonor];
                //LOG << " " << j << "-th nearest neighbor is donor cell " << iDonor << " " << donorPoly;
                std::list<Polygon_with_holes_2> pwhList;
                CGAL::intersection(selfPoly, donorPoly, std::back_inserter(pwhList));
                if (!pwhList.empty())
                {
                    double area = 0.0;
                    for (auto pwh : pwhList)
                    {
                        area += std::abs(CGAL::to_double(pwh.outer_boundary().area()));
                    }
                    //LOG << " intersects donor cell " << j << "(" << mDonorCellFaceAreas[j] << "), union = " << area;
                    mCellFaceIntersections[i].push_back(CellFaceIntersection(area, iDonor));
                    areaIntersected += area;
                    overlap = areaIntersected / selfArea;
                    if (overlap >= overlapMin)
                        break;
                }
                //LOG << std::endl;
            }
            //LOG << " overlap " << overlap * 100.0 << "%";
            if (overlap < overlapMin)
            {
                //LOG << " **** INSUFFICIENT OVERLAP ****" << std::endl;
            }
            else
            {
                //LOG << std::endl;
            }

            kStart = k;
        }
    }

    annDeallocPt(selfCentroid);
    delete[] nnIdx;
    delete[] dists;
}

void
SlidingInterface::MapMesh(const IterationContext& iteration)
{
    // FIXME: the latter part of InitializeMapping should go here
}

void
SlidingInterface::MapData(const Model& model, InterfaceDataAdaptorBase* adaptor) const
{
    std::ostream& LOG = Communicator::GetInstance()->Console();
#if SLIDING_INTERFACE_DEBUG
    LOG << "MapData:" << std::endl;
#endif

    for (size_t i = 0; i < mSelfCellFaceIndices.size(); ++i)
    {
        const PatchCellFaceIndex& pcfi = mSelfCellFaceIndices[i];
        const BlockPatch& bp = Roster::GetInstance()->GetBlockPatch(pcfi.BlockPatchID);
        Structured<double>& U = adaptor->GetBlockData(bp.BlockID());
        size_t dof = U.DOF();

        IndexIJK i3ghost = bp.I3Ghost();
        IndexIJK ijkCellFace(pcfi.I, pcfi.J, pcfi.K);
        IndexIJK ijk = ijkCellFace + i3ghost;

        double area = mSelfCellFaceAreas[i];
        double* u = U(ijk);
        for (size_t l = 0; l < dof; ++l)
            u[l] = 0.0;
#if SLIDING_INTERFACE_DEBUG
        LOG << "Cell " << i << ", " << bp.BlockID() << ':' << ijk << " CellFace " << ijkCellFace << std::endl;
#endif

        const std::list<CellFaceIntersection>& cfints = mCellFaceIntersections[i];
        for (const auto& cfint : cfints)
        {
            const PatchCellFaceIndex& pcfiD = mDonorCellFaceIndices[cfint.CellFaceIndex];
            const BlockPatch& bpD = Roster::GetInstance()->GetBlockPatch(pcfiD.BlockPatchID);
            const Structured<double>& UD = adaptor->GetBlockPatchData(pcfiD.BlockPatchID);
            double areaInt = cfint.Area;
            IndexIJK ijkCellFaceD = IndexIJK(pcfiD.I, pcfiD.J, pcfiD.K);
            IndexIJK ijkD = ijkCellFaceD + bpD.I3Interior();

#if SLIDING_INTERFACE_DEBUG
            LOG << " Donor " << bpD.BlockID() << ':' << ijkD << " CellFace " << ijkCellFaceD << ", areaInt " << areaInt << std::endl;
#endif
            double* ud = UD(ijkD);
            for (size_t l = 0; l < dof; ++l)
                u[l] += ud[l] * areaInt / area;
        }
    }

#if SLIDING_INTERFACE_DEBUG
    LOG << "MapData: data copied, before conversion to local frame" << std::endl;
    DumpInterfaceGhostCells(LOG, model);
#endif

    for (const auto& bp : SelfBlockPatches())
    {
        Structured<double>& U = adaptor->GetBlockData(bp.BlockID());
        ConvertMappedDataToLocalFrame(model, U, bp);
    }

#if SLIDING_INTERFACE_DEBUG
    LOG << "MapData: data copied, after conversion to local frame" << std::endl;
    DumpInterfaceGhostCells(LOG, model);
#endif
}

#include <fstream>

void write_polygon(std::ostream& o, const SlidingInterface::Polygon_2& poly)
{
    for (SlidingInterface::Polygon_2::Vertex_const_iterator i = poly.vertices_begin();
        i != poly.vertices_end(); ++i)
    {
        SlidingInterface::Point_2 pt = *i;
        o << CGAL::to_double(pt.x()) << '\t' << CGAL::to_double(pt.y()) << std::endl;
    }
    o << CGAL::to_double(poly[0].x()) << '\t' << CGAL::to_double(poly[0].y()) << std::endl;
}

void
SlidingInterface::DumpIntersection(size_t iSelf, ANNidxArray nnIdx, int k) const
{
    std::ofstream dump("overlap_failed.dat");
    write_polygon(dump, mSelfCellFaces[iSelf]);
    dump << std::endl << std::endl;
    for (size_t i = 0; i < k; ++i)
    {
        write_polygon(dump, mDonorCellFaces[nnIdx[i]]);
        dump << std::endl;
    }
    dump.close();

    dump.open("overlap_failed_alldonors.dat");
    write_polygon(dump, mSelfCellFaces[iSelf]);
    dump << std::endl << std::endl;
    for (size_t i = 0; i < mDonorCellFaces.size(); ++i)
    {
        write_polygon(dump, mDonorCellFaces[i]);
        dump << std::endl;
    }
    dump.close();
}

