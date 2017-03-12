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

#include "Roster.h"
#include "PlanarMapping.h"
#include "IndexUtils.h"
#include "Vector3.h"
#include "Communicator.h"

PlanarMapping::PlanarMapping(const Structured<double>& patchXYZ, const BlockPatch& bp)
:   mXYZ(patchXYZ), mBP(bp)
{
    IndexRange cfr = mBP.CellFaceRange();
    mCells.Allocate(4, cfr);
    mCells = -1;
    mD.Allocate(1, cfr);
    mD = -1.0;
    ComputeCentroids(mBP, mC, mXYZ);
}

PlanarMapping::~PlanarMapping()
{
    delete[] mCells.Data;
    delete[] mD.Data;
}

bool Contains(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p)
{
    // (3)
    //  | \    .
    //  |  \   .
    //  |   \  .
    // (1)--(2)
    Vector3 v12 = p2 - p1;
    Vector3 v13 = p3 - p1;
    Vector3 v1p = p - p1;
    double u = dot_product(v1p, v12) / v12.MagSq();
    double v = dot_product(v1p, v13) / v13.MagSq();
//std::cout << "TriContains: v12 = " << v12 << ", v13 = " << v13 << ", v1p = " << v1p << std::endl;
//std::cout << "TriContains: u = " << u << ", v = " << v << ", u + v = " << u + v << std::endl;
    return 0.0 <= u && u <= 1.0 && 0.0 <= v && v <= 1.0 && u + v <= 1.0;
}

bool Contains(const Structured<double>& xyz, const IndexIJK& i, const IndexIJK& i1, const IndexIJK& i2, const Vector3& p, IndexIJK& icell, bool debug = false)
{
    //std::cout << "Contains: " << i << ", " << i1 << ", " << i2 << ", " << p << std::endl;
    //  (NW)--(N)--(NE)
    //    |    |    |
    //   (W)--(C)--(E)
    //    |    |    |
    //  (SW)--(S)--(SE)
    //IndexRange r = xyz.GetRange();

    Vector3 pC(xyz(i));
    Vector3 pE(xyz(i + i1));
    Vector3 pW(xyz(i - i1));
    Vector3 pN(xyz(i + i2));
    Vector3 pS(xyz(i - i2));
    Vector3 pNE(xyz(i + i1 + i2));
    Vector3 pNW(xyz(i - i1 + i2));
    Vector3 pSE(xyz(i + i1 - i2));
    Vector3 pSW(xyz(i - i1 - i2));

    bool CSE, CEN, CNW, CWS, SSEE, ENEN, WNNW, WSWS;
    CSE = Contains(pC, pS, pE, p);
    CEN = Contains(pC, pE, pN, p);
    CNW = Contains(pC, pN, pW, p);
    CWS = Contains(pC, pW, pS, p);
    SSEE = Contains(pS, pSE, pE, p);
    ENEN = Contains(pE, pNE, pN, p);
    WNNW = Contains(pW, pN, pNW, p);
    WSWS = Contains(pW, pSW, pS, p);

    bool contains = CSE || CEN || CNW || CWS || SSEE || ENEN || WNNW || WSWS;

    if (CSE || SSEE) icell = IndexIJK(i + i1);
    else if (CEN || ENEN) icell = IndexIJK(i + i1 + i2);
    else if (CNW || WNNW) icell = IndexIJK(i + i2);
    else if (CWS || WSWS) icell = IndexIJK(i);

    if (contains && debug)
    {
        std::cout << "Contains: " << i << ", " << i1 << ", " << i2 << ", " << p << std::endl;
        std::cout << i << ":" << CEN << ENEN << CNW << WNNW << CWS << WSWS << CSE << SSEE << std::endl;
        std::cout << p << std::endl << pC << std::endl << pE << std::endl << pW << std::endl << pN << std::endl << pS << std::endl << pNE << std::endl << pNW << std::endl << pSE << std::endl << pSW << std::endl;
        assert(false);
    }

    return contains;
}

bool Contains(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4, const Vector3& p)
{
    // (p3)---(p4)
    //  |       |
    // (p1)---(p2)

    bool t1, t2;
    t1 = Contains(p1, p2, p3, p);
    t2 = Contains(p2, p4, p3, p);
    return t1 || t2;
}

void
PlanarMapping::ComputeCentroids(const BlockPatch& bp, Structured<double>& centroids, const Structured<double>& XYZ) const
{
    IndexRange cfr = bp.CellFaceRange();
    IndexIJK i1 = bp.I1();
    IndexIJK i2 = bp.I2();

    if (centroids.Data != NULL)
    {
        delete[] centroids.Data;
    }
    centroids.Allocate(4, cfr); // centroid position + cell size metric (currently diagonal length)

    for (IndexIterator i(cfr); !i.IsEnd(); i.Advance())
    {
        IndexIJK ijk = i.Index();
        Vector3 p1(XYZ(ijk - i1 - i2));
        Vector3 p2(XYZ(ijk - i2));
        Vector3 p3(XYZ(ijk - i1));
        Vector3 p4(XYZ(ijk));
        Vector3 p = 0.25 * (p1 + p2 + p3 + p4);
        double* c = centroids(ijk);
        c[0] = p[0];
        c[1] = p[1];
        c[2] = p[2];

        // take the longer diagonal length (***squared***) to be the cell size metric
        double d1, d2;
        d1 = (p4 - p1).MagSq();
        d2 = (p3 - p2).MagSq();
        c[3] = std::max(d1, d2);
    }
}

void
PlanarMapping::InterpolateOn(double& c12, double& c13, const Vector3& p, const Structured<double>& XYZ,
    const IndexIJK& ijk, const IndexIJK& i1, const IndexIJK& i2) const
{
    // 3---4    ijk points to the cell center index "o"
    // | o |    i1 is the first index basis for this plane.
    // 1---2    i2 is the second index basis for this plane.

    Vector3 p1(XYZ(ijk - i1 - i2));
    Vector3 p2(XYZ(ijk - i2));
    Vector3 p3(XYZ(ijk - i1));
    Vector3 p4(XYZ(ijk));
    Vector3 p0(p1, p);

    Vector3 e12(p1, p2);
    Vector3 e13(p1, p3);

    c12 = dot_product(p0, e12) / e12.MagSq();
    c13 = dot_product(p0, e13) / e13.MagSq();
}

void
PlanarMapping::InitializeMapping()
{
    mD = -1.0;
}

#if 1
void
PlanarMapping::MapOn(const BlockPatch& bpd, const Structured<double>& XYZDonor, const Vector3& angleSelf, const Vector3& angleDonor, bool debug)
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    IndexRange cfr, cfrd;
    IndexIJK i1, i2, i1d, i2d;
    cfr = mBP.CellFaceRange();
    i1 = mBP.I1();
    i2 = mBP.I2();
    cfrd = bpd.CellFaceRange();
    i1d = bpd.I1();
    i2d = bpd.I2();

    // donor mesh centroids
    Structured<double> centroidsD;
    ComputeCentroids(bpd, centroidsD, XYZDonor);

    double RSelf[3][3], RDonor[3][3];
    Vector3::RotationMatrix(RSelf, angleSelf);
    Vector3::RotationMatrix(RDonor, angleDonor);

    for (IndexIterator i(cfr); !i.IsEnd(); i.Advance())
    {
        IndexIJK ijk = i.Index();
        Vector3 centroid(mC(ijk));
        centroid.Apply(RSelf);
        double cellLengthSq = mC(ijk)[3];
        double* distSq = mD(ijk);
        for (IndexIterator j(cfrd); !j.IsEnd(); j.Advance())
        {
            IndexIJK ijkd = j.Index();
            Vector3 centroidD(centroidsD(ijkd));
            centroidD.Apply(RDonor);
            double dsq = (centroid - centroidD).MagSq();
            //if (dsq < cellLengthSq * 1.5 && (*distSq < 0.0 || dsq < *distSq))
            if (*distSq < 0.0 || dsq < *distSq)
            {
                double c12, c13;
                // FIXME: the next line, cannot be right, coordinates in XYZDonor are not rotated properly.
                // FIXME: but after all we are not using its results (c12 and c13), look at the if statement that follows.
                InterpolateOn(c12, c13, centroid, XYZDonor, ijkd, i1d, i2d);
                if (true || std::abs(c12) <= 1.0 && std::abs(c13) <= 1.0)
                {
                    *distSq = dsq;
                    int* cell = mCells(ijk);
                    cell[0] = bpd.BlockID();
                    cell[1] = ijkd.I;
                    cell[2] = ijkd.J;
                    cell[3] = ijkd.K;
                    if (debug)
                    {
                        LOG << "PlanarMapping: " << ijk << " -> Block " << cell[0] << ": " << ijkd << std::endl;
                    }
                }
            }
        }
    }

    delete[] centroidsD.Data;
}
#else
void
PlanarMapping::MapOn(int blockID, const Structured<double>& XYZDonor, const Vector3& angle)
{
    Vector3 e1, e2;
    PatchBasis(e1, e2, mXYZ);
    //std::cout << "Patch basis - " << e1 << ", " << e2 << std::endl;

    IndexRange mr = mXYZ.GetRange();
    //std::cout << "MeshRange = " << mr << std::endl;
    IndexRange cfr;
    IndexIJK i1, i2, i3;
    CellFaceRange(cfr, i1, i2, i3, mr);
    IndexRange mrDonor = XYZDonor.GetRange();
    //std::cout << "DonorMeshRange = " << mrDonor << std::endl;
    IndexRange cfrDonor;
    IndexIJK i1d, i2d, i3d;
    CellFaceRange(cfrDonor, i1d, i2d, i3d, mrDonor);

    double R[3][3];
    Vector3::RotationMatrix(R, angle);

    for (IndexIterator i(cfr); !i.IsEnd(); i.Advance())
    {
        IndexIJK ijk = i.Index();
        Vector3 p1(mXYZ(ijk + i3 - i1 - i2));
        Vector3 p2(mXYZ(ijk + i3 - i2));
        Vector3 p3(mXYZ(ijk + i3 - i1));
        Vector3 p4(mXYZ(ijk + i3));
        Vector3 p = 0.25 * (p1 + p2 + p3 + p4);
        //double* D = mD(ijk);
        int* cell = mCells(ijk);
        for (IndexIterator j(cfrDonor); !j.IsEnd(); j.Advance())
        {
            IndexIJK ijkd = j.Index();
            Vector3 p1d(XYZDonor(ijkd + i3d - i1d - i2d));
            Vector3 p2d(XYZDonor(ijkd + i3d - i2d));
            Vector3 p3d(XYZDonor(ijkd + i3d - i1d));
            Vector3 p4d(XYZDonor(ijkd + i3d));
            p1d = p1d.Apply(R);
            p2d = p2d.Apply(R);
            p3d = p3d.Apply(R);
            p4d = p4d.Apply(R);
            bool contains = Contains(p1d, p2d, p3d, p4d, p);
            if (contains)
            {
                //std::cout << ijk << " -> " << ijkd << ", " << p1d << p2d << p3d << p4d << p << std::endl;
                cell[0] = blockID;
                cell[1] = ijkd[0];
                cell[2] = ijkd[1];
                cell[3] = ijkd[2];
                break;
            }
        }
    }
}
#endif

bool
PlanarMapping::AllMapped() const
{
    bool allMapped = true;
    for (IndexIterator i(mCells.GetRange()); !i.IsEnd(); i.Advance())
    {
        int* cell = mCells(i.Index());
        if (cell[0] < 0)
        {
            allMapped = false;
            break;
        }
    }
    return allMapped;
}

Vector3
PlanarMapping::PatchNormal(const Structured<double>& XYZ) const
{
    Vector3 e1, e2;
    PatchBasis(e1, e2, XYZ);
    Vector3 normal = cross_product(e1, e2);
    normal.Normalize();
    return normal;
}

void
PlanarMapping::PatchBasis(Vector3& e1, Vector3& e2, const Structured<double>& XYZ) const
{
    // Find the normal vector of the interface plane.
    IndexRange bmr = Roster::GetInstance()->GetBlock(mBP.BlockID())->MeshRange();
    IndexRange mr = XYZ.GetRange();
    Vector3 p1(XYZ(mr.Start)), p2(XYZ(mr.End));
    IndexIJK i3(mr.Start), i4(mr.End);
    Direction dir = IndexUtils::PatchDirection(bmr, mr);
    switch (dir)
    {
    case I:
    case INEG:
        i3.J = mr.End.J;
        i4.J = mr.Start.J;
        break;
    case J:
    case JNEG:
        i3.I = mr.End.I;
        i4.I = mr.Start.I;
        break;
    case K:
    case KNEG:
        i3.J = mr.End.J;
        i4.J = mr.Start.J;
        break;
    }
    Vector3 p3(XYZ(i3)), p4(XYZ(i4));
    Vector3 v1 = p2 - p1, v2 = p4 - p3;
    // Gram-Schmidt orthogonalization
    e1 = v1.Normalized();
    e2 = v2 - dot_product(v1, v2) * e1;
    e2.Normalize();
}

void
PlanarMapping::MapData(const Model& model, Structured<double>& U, const BlockPatch& bpd, const Structured<double>& UDonor) const
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    int ndof = U.DOF();
    assert(ndof == UDonor.DOF());

    IndexIJK i3ghost = mBP.I3Ghost();
    IndexIJK i3interiorD = bpd.I3Interior();

    for (IndexIterator it = mCells.GetRange(); !it.IsEnd(); it.Advance())
    {
        IndexIJK ijk = it.Index();
        int* cell = mCells(ijk);
        if (cell[0] == bpd.BlockID())
        {
            IndexIJK ijkd(cell[1], cell[2], cell[3]);
            IndexIJK ijkSelf = ijk + i3ghost, ijkDonor = ijkd + i3interiorD;
            double* u = U(ijkSelf);
            double* ud = UDonor(ijkDonor);
            for (int l = 0; l < ndof; ++l)
            {
                u[l] = ud[l];
            }
#if 0
            log << "Mapped " << ijkSelf << "(self) = " << ijkDonor << "(donor):";
            for (int l = 0; l < ndof; ++l)
                log << ' ' << u[l];
            log << " : ";
            for (int l = 0; l < ndof; ++l)
                log << ' ' << ud[l];
            log << std::endl;
#endif
        }
    }
}

