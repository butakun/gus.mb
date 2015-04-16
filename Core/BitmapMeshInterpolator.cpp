#include "BitmapMeshInterpolator.h"
#include "Matrix33.h"
#include "GeometryTools.h"
#undef NDEBUG
#include <cassert>

void
CellLocator::AddPatch(const BlockPatch& bp, const Structured<double>& XYZ)
{
    assert(mPatchMeshMap.find(bp.UniqueID()) == mPatchMeshMap.end());

    mBlockPatches.push_back(bp);
    mPatchMeshMap[bp.UniqueID()] = XYZ;
}

bool
CellLocator::LocatePoint(PatchCell& cell, const Vector3& p) const
{
    for (BlockPatches::const_iterator i = mBlockPatches.begin();
        i != mBlockPatches.end(); ++i)
    {
        const BlockPatch& bp = *i;
        const Structured<double>& XYZ = mPatchMeshMap.at(bp.UniqueID());
        IndexIJK ijk;
        if (LocatePointInPatch(ijk, p, bp, XYZ))
        {
            cell = PatchCell(bp.UniqueID(), ijk);
            return true;
        }
    }
    return false;
}

bool
CellLocator::LocatePointInPatch(IndexIJK& ijk, const Vector3& p, const BlockPatch& bp, const Structured<double>& XYZ) const
{
    IndexRange cfr = bp.CellFaceRange();
    IndexIJK i1 = bp.I1();
    IndexIJK i2 = bp.I2();

    for (IndexIterator i(cfr); !i.IsEnd(); i.Advance())
    {
        IndexIJK ijk_ = i.Index();
        Vector3 p1(XYZ(ijk_ - i1 - i2));
        Vector3 p2(XYZ(ijk_ - i2));
        Vector3 p3(XYZ(ijk_));
        Vector3 p4(XYZ(ijk_ - i1));
        bool inside = IsPointInside(p, p1, p2, p3, p4);
        if (inside)
        {
            ijk = ijk_;
            return inside;
        }
    }
    return false;
}

bool
CellLocator::IsPointInside(const Vector3& p, const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4) const
{
    // p4 --- p3
    // |   p   |
    // p1 --- p2
    return IsPointInsideTriangle(p, p1, p2, p3) || IsPointInsideTriangle(p, p1, p3, p4);
}

bool
CellLocator::IsPointInsideTriangle(const Vector3& p, const Vector3& p1, const Vector3& p2, const Vector3& p3) const
{
    //    p3
    //   /p|
    // p1--p2
    Vector3 n = cross_product(p2 - p1, p3 - p1);
    n.Normalize();
    Vector3 pp1 = p - p1;
    Vector3 pp = p1 + pp1 - dot_product(pp1, n) * n;

    Matrix33 M;
    M(0, 0) = p1.X();
    M(1, 0) = p1.Y();
    M(2, 0) = p1.Z();
    M(0, 1) = p2.X();
    M(1, 1) = p2.Y();
    M(2, 1) = p2.Z();
    M(0, 2) = p3.X();
    M(1, 2) = p3.Y();
    M(2, 2) = p3.Z();

    Matrix33 Minv = M.Inverse();
    Vector3 A = Minv * pp;

#if 1
    return A[0] >= 0.0 && A[1] >= 0.0 && A[2] >= 0.0 && (A[0] + A[1] + A[2]) <= 1.00001;
#else
    bool result = A[0] >= 0.0 && A[1] >= 0.0 && A[2] >= 0.0 && (A[0] + A[1] + A[2]) <= 1.00001;
    std::cout << "IsPointInsideTriangle: p = " << p << ", p1p2p3 = " << p1 << " - " << p2 << " - " << p3 << " --> " << result << std::endl;
    std::cout << "pp = " << pp << ", A = " << A[0] << " + " << A[1] << " + " << A[2] << " = " << A[0] + A[1] + A[2] << std::endl;
    std::cout << pp.Y() << '\t' << pp.Z() << std::endl << std::endl << std::endl;
    std::cout << p1.Y() << '\t' << p1.Z() << std::endl;
    std::cout << p2.Y() << '\t' << p2.Z() << std::endl;
    std::cout << p3.Y() << '\t' << p3.Z() << std::endl;
    return result;
#endif
}

BitmapMeshInterpolator::BitmapMeshInterpolator()
{
}

BitmapMeshInterpolator::~BitmapMeshInterpolator()
{
}

void
BitmapMeshInterpolator::AddPatch(int which, const BlockPatch& bp, const Structured<double>& XYZ)
{
    assert(which == 0 || which == 1);

    if (which == 0)
    {
        mPatches1.push_back(bp);
    }
    else
    {
        mPatches2.push_back(bp);
    }
    mPatchMeshes[bp.UniqueID()] = XYZ;
}

void
BitmapMeshInterpolator::GenerateBitmapPixels()
{
    // Count the number of cells in each group of patches
    int cellCount1 = 0;
    for (BlockPatches::const_iterator i = mPatches1.begin();
        i != mPatches1.end(); ++i)
    {
        IndexIJK shape = i->CellRange().Shape();
        cellCount1 += (shape.I + shape.J + shape.K);
    }

    int cellCount2 = 0;
    for (BlockPatches::const_iterator i = mPatches2.begin();
        i != mPatches2.end(); ++i)
    {
        IndexIJK shape = i->CellRange().Shape();
        cellCount2 += (shape.I + shape.J + shape.K);
    }

    // We choose the denser group as the reference for pixel generation
    size_t refPatchIndex = 0; // 0 if refpatch -> patches1, 1 otherwise
    size_t donorPatchIndex = 1;
    BlockPatches* refPatches = NULL;
    BlockPatches* donorPatches = NULL;
    if (cellCount2 <= cellCount1)
    {
        refPatchIndex = 0;
        refPatches = &mPatches1;
        donorPatchIndex = 1;
        donorPatches = &mPatches2;
    }
    else
    {
        refPatchIndex = 1;
        refPatches = &mPatches2;
        donorPatchIndex = 0;
        donorPatches = &mPatches1;
    }

    // Generate pixels
    //int NSUBDIV = 10;
    int NSUBDIV = 5;
    Structured<double> XYZSub(3, IndexRange(IndexIJK(0, 0, 0), IndexIJK(NSUBDIV, NSUBDIV, 0)));

    double patchArea = 0.0;
    for (BlockPatches::const_iterator i = refPatches->begin();
        i != refPatches->end(); ++i)
    {
        const BlockPatch& bp = *i;
        IndexIJK i1 = bp.I1();
        IndexIJK i2 = bp.I2();
        assert(mPatchMeshes.find(bp.UniqueID()) != mPatchMeshes.end());
        Structured<double>& XYZ = mPatchMeshes.at(bp.UniqueID());
        for (IndexIterator itor(bp.CellFaceRange()); !itor.IsEnd(); itor.Advance())
        {
            IndexIJK icell = itor.Index();
            // 4 --- 3
            // |     |
            // 1 --- 2
            Vector3 p1, p2, p3, p4;
            p1 = Vector3(XYZ(icell - i1 - i2));
            p2 = Vector3(XYZ(icell      - i2));
            p3 = Vector3(XYZ(icell          ));
            p4 = Vector3(XYZ(icell - i1     ));
            double cellArea = GeometryTools::QuadArea(p1, p2, p3, p4);
            patchArea += cellArea;
            for (int is = 0; is <= NSUBDIV; ++is)
            {
                double xi = double(is) / double(NSUBDIV);
                for (int js = 0; js <= NSUBDIV; ++js)
                {
                    double eta = double(js) / double(NSUBDIV);
                    Vector3 p = (1.0 - xi) * (1.0 - eta) * p1 + xi * (1.0 - eta) * p2 + (1.0 - xi) * eta * p4 + xi * eta * p3;
                    double* v = XYZSub(is, js, 0);
                    v[0] = p.X();
                    v[1] = p.Y();
                    v[2] = p.Z();
                }
            }
            for (int is = 0; is < NSUBDIV; ++is)
            {
                for (int js = 0; js < NSUBDIV; ++js)
                {
                    Vector3 p1, p2, p3, p4;
                    p1 = Vector3(XYZSub(is    , js    , 0));
                    p2 = Vector3(XYZSub(is + 1, js    , 0));
                    p3 = Vector3(XYZSub(is + 1, js + 1, 0));
                    p4 = Vector3(XYZSub(is    , js + 1, 0));
                    Vector3 pixel = 0.25 * (p1 + p2 + p3 + p4);
                    mPixels.push_back(pixel);
                    mPixelWeights.push_back(cellArea / double(NSUBDIV * NSUBDIV));
                    PatchIndex pi;
                    pi.PatchID[refPatchIndex] = bp.UniqueID();
                    pi.PatchIJK[refPatchIndex] = icell;
                    mPixelMapping.push_back(pi);
                }
            }
        }
    }
    delete[] XYZSub.Data;

    // Scale the weight by the total area of the patches
    for (std::vector<double>::iterator i = mPixelWeights.begin();
        i != mPixelWeights.end(); ++i)
    {
        *i /= patchArea;
    }

    // Match the pixels onto the donor patches
    CellLocator cellLocator;
    for (BlockPatches::const_iterator i = donorPatches->begin();
        i != donorPatches->end(); ++i)
    {
        const BlockPatch& bp = *i;
        assert(mPatchMeshes.find(bp.UniqueID()) != mPatchMeshes.end());
        Structured<double>& XYZ = mPatchMeshes.at(bp.UniqueID());
        cellLocator.AddPatch(bp, XYZ);
    }

    assert(mPixels.size() == mPixelMapping.size());
    for (size_t i = 0; i < mPixels.size(); ++i)
    {
        const Vector3& pixel = mPixels[i];
        CellLocator::PatchCell patchCell;
        bool found = cellLocator.LocatePoint(patchCell, pixel);
        assert(found);
        mPixelMapping[i].PatchID[donorPatchIndex] = patchCell.id;
        mPixelMapping[i].PatchIJK[donorPatchIndex] = patchCell.index;
    }
}

void
BitmapMeshInterpolator::ComputeCellFaceWeights()
{
    assert(mCellFaceWeights.empty());

    // GenerateBitmapPixels has bee completed.
    for (int ibps = 0; ibps < 2; ++ibps)
    {
        BlockPatches* bps = ibps == 0 ? &mPatches1 : &mPatches2;
        for (BlockPatches::const_iterator ibp = bps->begin();
            ibp != bps->end(); ++ibp)
        {
            const BlockPatch& bp = *ibp;
            IndexRange cfr = bp.CellFaceRange();
            Structured<double> cfw(1, cfr);
            cfw = 0.0;
            mCellFaceWeights[bp.UniqueID()] = cfw;
            Structured<size_t> cfm(1, cfr);
            mCellFaceMapping[bp.UniqueID()] = cfm;
        }
    }

    for (size_t i = 0; i < mPixels.size(); ++i)
    {
        double pixelWeight = mPixelWeights[i];
        const PatchIndex& pi = mPixelMapping[i];
        for (int ip = 0; ip < 2; ++ip)
        {
            assert(mCellFaceWeights.find(pi.PatchID[ip]) != mCellFaceWeights.end());
            Structured<double>& cfw = mCellFaceWeights[pi.PatchID[ip]];
            cfw(pi.PatchIJK[ip])[0] += pixelWeight;
        }
    }
}

