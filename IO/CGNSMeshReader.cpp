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
// $Id: CGNSMeshReader.cpp 321 2014-08-26 09:52:44Z kato $

#include "CGNSMeshReader.h"
#include "Block.h"
#include "IndexUtils.h"
namespace CGNS {
#include "cgnslib.h"
}
#include <iostream>
#include <cassert>

const size_t CGNS_NAME_LENGTH = 33;

typedef struct {
    CGNS::CG_BCType_t type;
    const char* name;
} bocotype_to_string_t;

bocotype_to_string_t bocotype_to_string[] = {
    { (CGNS::CG_BCType_t)CG_UserDefined, "UserDefined" },
    { CGNS::CG_BCAxisymmetricWedge, "AxisymmetricWedge" },
    { CGNS::CG_BCDegenerateLine, "DegenerateLine" },
    { CGNS::CG_BCExtrapolate, "Extrapolate" },
    { CGNS::CG_BCDegeneratePoint, "DegeneratePoint" },
    { CGNS::CG_BCDirichlet, "Dirichlet" },
    { CGNS::CG_BCFarfield, "Farfield" },
    { CGNS::CG_BCNeumann, "Neumann" },
    { CGNS::CG_BCGeneral, "General" },
    { CGNS::CG_BCInflow, "Inflow" },
    { CGNS::CG_BCOutflow, "Outflow" },
    { CGNS::CG_BCInflowSubsonic, "InflowSubsonic" },
    { CGNS::CG_BCOutflowSubsonic, "OutflowSubsonic" },
    { CGNS::CG_BCInflowSupersonic, "InflowSupersonic" },
    { CGNS::CG_BCOutflowSupersonic, "OutflowSupersonic" },
    { CGNS::CG_BCSymmetryPlane, "SymmetryPlane" },
    { CGNS::CG_BCTunnelInflow, "TunnelInflow" },
    { CGNS::CG_BCSymmetryPolar, "SymmetryPolar" },
    { CGNS::CG_BCTunnelOutflow, "TunnelOutflow" },
    { CGNS::CG_BCWallViscous, "WallViscous" },
    { CGNS::CG_BCWall, "Wall" },
    { CGNS::CG_BCWallViscousHeatFlux, "WallViscousHeatFlux" },
    { CGNS::CG_BCWallInviscid, "WallInviscid" },
    { CGNS::CG_BCWallViscousIsothermal, "WallViscousIsothermal" },
    { CGNS::CG_FamilySpecified, "FamilySpecified" },
    { (CGNS::CG_BCType_t)0, "Null" }
};

static void cgns_assert(int err)
{
    if (err == 0)
        return;
    std::cerr << CGNS::cg_get_error() << std::endl;
    assert(false);
}

const char* bocotype_string(CGNS::CG_BCType_t type)
{
    bocotype_to_string_t* p = bocotype_to_string;
    for (; p->type != 0; ++p)
        if (p->type == type)
            return p->name;
    return "Null";
}

void canonize_transform(IndexRange& selfRange, IndexRange& donorRange, int* transform)
{
    if (selfRange.IsCanonical())
    {
        return;
    }

    if (selfRange.Start.I > selfRange.End.I)
    {
        int i = std::abs(transform[0]) - 1;
        int tmp;
        tmp = donorRange.Start[i];
        donorRange.Start[i] = donorRange.End[i];
        donorRange.End[i] = tmp;
    }
    if (selfRange.Start.J > selfRange.End.J)
    {
        int i = std::abs(transform[1]) - 1;
        int tmp;
        tmp = donorRange.Start[i];
        donorRange.Start[i] = donorRange.End[i];
        donorRange.End[i] = tmp;
    }
    if (selfRange.Start.K > selfRange.End.K)
    {
        int i = std::abs(transform[2]) - 1;
        int tmp;
        tmp = donorRange.Start[i];
        donorRange.Start[i] = donorRange.End[i];
        donorRange.End[i] = tmp;
    }

    selfRange.Canonize();
}

CGNSMeshReader::CGNSMeshReader(const char* filename)
:   mFileName(filename)
{
    ReadBlockInfo();
}

CGNSMeshReader::~CGNSMeshReader()
{
}

void
CGNSMeshReader::ReadBlockInfo()
{
    int err;

    int fn;
    err = CGNS::cg_open(mFileName.c_str(), CG_MODE_READ, &fn);
    cgns_assert(err);

    int nbases;
    err = CGNS::cg_nbases(fn, &nbases);
    cgns_assert(err);
    assert(nbases > 0);
    if (nbases > 1)
    {
        std::cerr << "Warning: there are " << nbases << " bases in the CGNS file " << mFileName
            << " but we assume the first base contains the mesh." << std::endl;
    }

    int B = 1;
    char basename[CGNS_NAME_LENGTH];
    int cell_dim, phys_dim;
    err = CGNS::cg_base_read(fn, B, basename, &cell_dim, &phys_dim);
    cgns_assert(err);

    std::cout << "cell_dim = " << cell_dim << ", phys_dim = " << phys_dim << std::endl;
    assert(cell_dim == 3);
    assert(phys_dim == 3);

    int nzones;
    err = CGNS::cg_nzones(fn, B, &nzones);
    cgns_assert(err);

    std::cout << "There are " << nzones << " zones." << std::endl;

    mDomain = DomainInfo(); // initialize with an empty domain info.
    for (int Z = 1; Z <= nzones; ++Z)
    {
        CGNS::CG_ZoneType_t zoneType;
        err = CGNS::cg_zone_type(fn, B, Z, &zoneType);
        cgns_assert(err);
        assert(zoneType == CGNS::CG_Structured);

        char zoneName[CGNS_NAME_LENGTH];
        CGNS::cgsize_t size[9];
        err = CGNS::cg_zone_read(fn, B, Z, zoneName, size);
        cgns_assert(err);

        std::cout << "Zone " << Z << " Name " << zoneName << std::endl;
        std::cout << "NVertices: " << size[0] << "x" << size[1] << "x" << size[2] << std::endl;
        std::cout << "NCells: " << size[3] << "x" << size[4] << "x" << size[5] << std::endl;
        int imin = 0, jmin = 0, kmin = 0, imax = size[0] - 1, jmax = size[1] - 1, kmax = size[2] - 1;
        IndexRange meshRange(imin, jmin, kmin, imax, jmax, kmax);

        int ncoords;
        err = CGNS::cg_ncoords(fn, B, Z, &ncoords);
        cgns_assert(err);
        assert(ncoords == 3);

        BlockInfo bi(zoneName, Z, meshRange);
        mDomain.AddBlockInfo(bi);
    }

    for (int Z = 1; Z <= nzones; ++Z)
    {
        std::cout << "DEBUG: Z = " << Z << std::endl;
        BlockInfo& bi = mDomain.FindBlockInfo(Z);

        int nbocos;
        err = CGNS::cg_nbocos(fn, B, Z, &nbocos);
        cgns_assert(err);

        bi.GetBCInfos().clear();
        for (int BC = 1; BC <= nbocos; ++BC)
        {
            char boconame[CGNS_NAME_LENGTH];
            CGNS::CG_BCType_t bocotype;
            CGNS::CG_PointSetType_t ptset_type;
            CGNS::cgsize_t npnts;
            int NormalIndex[3];
            CGNS::cgsize_t NormalListSize;
            CGNS::CG_DataType_t NormalDataType;
            int ndataset;
            err = CGNS::cg_boco_info(fn, B, Z, BC, boconame, &bocotype, &ptset_type, &npnts, NormalIndex, &NormalListSize, &NormalDataType, &ndataset);
            cgns_assert(err);
            assert(ptset_type == CGNS::CG_PointRange);
            assert(npnts == 2);

            CGNS::cgsize_t pnts[2][3];
            err = CGNS::cg_boco_read(fn, B, Z, BC, &pnts[0][0], NULL);
            cgns_assert(err);

            IndexRange meshRange(
                pnts[0][0] - 1, pnts[0][1] - 1, pnts[0][2] - 1,
                pnts[1][0] - 1, pnts[1][1] - 1, pnts[1][2] - 1
                );

            bi.GetBCInfos().push_back(BCInfo(boconame, bocotype_string(bocotype), meshRange));
        }

        int n1to1, nconns;
        err = CGNS::cg_n1to1(fn, B, Z, &n1to1);
        cgns_assert(err);
        err = CGNS::cg_nconns(fn, B, Z, &nconns);
        cgns_assert(err);

        for (int I = 1; I <= n1to1; ++I)
        {
            char connName[CGNS_NAME_LENGTH];
            char donorName[CGNS_NAME_LENGTH];
            CGNS::cgsize_t range[6], donorRange[6];
            int transform[3];
            err = CGNS::cg_1to1_read(fn, B, Z, I, connName, donorName, range, donorRange, transform);
            cgns_assert(err);

            IndexRange selfMeshRange(
                range[0] - 1, range[1] - 1, range[2] - 1,
                range[3] - 1, range[4] - 1, range[5] - 1
                );
            int donorZone = mDomain.FindBlockInfo(donorName).Zone();
            IndexRange donorMeshRange(
                donorRange[0] - 1, donorRange[1] - 1, donorRange[2] - 1,
                donorRange[3] - 1, donorRange[4] - 1, donorRange[5] - 1
                );
            if (!selfMeshRange.IsCanonical())
            {
                canonize_transform(selfMeshRange, donorMeshRange, transform);
                std::cout << "Connectivity " << connName << " was not in canonical form, thus modified to be so." << std::endl;
            }
            Connectivity1to1Info conn(connName, Z, selfMeshRange, donorZone, donorMeshRange, transform);

            float rotCenter[3], rotAngle[3], trans[3];
            err = CGNS::cg_1to1_periodic_read(fn, B, Z, I, rotCenter, rotAngle, trans);
            if (err == CG_OK)
            {
                if (rotAngle[0] != 0.0 || rotAngle[1] != 0.0 || rotAngle[2] != 0.0)
                {
                    conn.SetPeriodicity(Connectivity1to1::ROTATION, Vector3(rotCenter), Vector3(rotAngle));
                }
            }

            bi.GetConn1to1s().push_back(conn);
        }
    }
}

void
CGNSMeshReader::ReadMesh(Structured<double>& XYZ, int Z) const
{
    //assert(XYZ.Data == NULL);

    int err;
    int fn, B = 1;
    GetZone(fn, B, Z);

    const BlockInfo& bi = mDomain.FindBlockInfo(Z);

    IndexRange mr = bi.MeshRange();
    int imin, jmin, kmin, imax, jmax, kmax;
    imin = mr.Start.I + 1;
    jmin = mr.Start.J + 1;
    kmin = mr.Start.K + 1;
    imax = mr.End.I + 1;
    jmax = mr.End.J + 1;
    kmax = mr.End.K + 1;

    CGNS::cgsize_t rangeMin[3] = {imin, jmin, kmin};
    CGNS::cgsize_t rangeMax[3] = {imax, jmax, kmax};

    // Temporary storage for the coord values
    double* xx = new double[imax * jmax * kmax];
    double* yy = new double[imax * jmax * kmax];
    double* zz = new double[imax * jmax * kmax];

    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateX", CGNS::CG_RealDouble, rangeMin, rangeMax, xx);
    cgns_assert(err);
    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateY", CGNS::CG_RealDouble, rangeMin, rangeMax, yy);
    cgns_assert(err);
    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateZ", CGNS::CG_RealDouble, rangeMin, rangeMax, zz);
    cgns_assert(err);

    Close(fn);

    // Allocates the storage. This will belong to the caller.
    if (XYZ.Data != NULL)
    {
        delete[] XYZ.Data;
    }
    XYZ.Data = new double[3 * imax * jmax * kmax];
    XYZ.SetRange(mr);
    XYZ.SetDOF(3);

    for (int k = 0; k < kmax; ++k)
    {
        for (int j = 0; j < jmax; ++j)
        {
            for (int i = 0; i < imax; ++i)
            {
                int ijk = k * jmax * imax + j * imax + i;
                XYZ(i, j, k)[0] = xx[ijk];
                XYZ(i, j, k)[1] = yy[ijk];
                XYZ(i, j, k)[2] = zz[ijk];
            }
        }
    }

    delete[] xx;
    delete[] yy;
    delete[] zz;
}

void
CGNSMeshReader::ReadPartialMesh(Structured<double>& XYZ, int Z, const IndexRange& range) const
{
    //assert(XYZ.Data == NULL);

    int err;
    int fn, B = 1;
    GetZone(fn, B, Z);

    const BlockInfo& bi = mDomain.FindBlockInfo(Z);

    IndexRange mr = bi.MeshRange();
    int imin, jmin, kmin, imax, jmax, kmax;
    imin = mr.Start.I + 1;
    jmin = mr.Start.J + 1;
    kmin = mr.Start.K + 1;
    imax = mr.End.I + 1;
    jmax = mr.End.J + 1;
    kmax = mr.End.K + 1;

    int iminRead, imaxRead, jminRead, jmaxRead, kminRead, kmaxRead;
    iminRead = range.Start.I + 1;
    jminRead = range.Start.J + 1;
    kminRead = range.Start.K + 1;
    imaxRead = range.End.I + 1;
    jmaxRead = range.End.J + 1;
    kmaxRead = range.End.K + 1;
    int idim, jdim, kdim;
    idim = imaxRead - iminRead + 1;
    jdim = jmaxRead - jminRead + 1;
    kdim = kmaxRead - kminRead + 1;

    CGNS::cgsize_t rangeMin[3] = {imin, jmin, kmin};
    CGNS::cgsize_t rangeMax[3] = {imax, jmax, kmax};

    // Temporary storage for the coord values
    double* xx = new double[imax * jmax * kmax];
    double* yy = new double[imax * jmax * kmax];
    double* zz = new double[imax * jmax * kmax];

    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateX", CGNS::CG_RealDouble, rangeMin, rangeMax, xx);
    cgns_assert(err);
    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateY", CGNS::CG_RealDouble, rangeMin, rangeMax, yy);
    cgns_assert(err);
    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateZ", CGNS::CG_RealDouble, rangeMin, rangeMax, zz);
    cgns_assert(err);

    Close(fn);

    // Allocates the storage. This will belong to the caller.
    if (XYZ.Data != NULL)
    {
        delete[] XYZ.Data;
    }
    XYZ.Allocate(3, IndexRange(0, 0, 0, idim - 1, jdim - 1, kdim - 1));

    for (int k = 0; k < kdim; ++k)
    {
        for (int j = 0; j < jdim; ++j)
        {
            for (int i = 0; i < idim; ++i)
            {
                int ijk = k * jmax * imax + j * imax + i;
                XYZ(i, j, k)[0] = xx[ijk];
                XYZ(i, j, k)[1] = yy[ijk];
                XYZ(i, j, k)[2] = zz[ijk];
            }
        }
    }

    delete[] xx;
    delete[] yy;
    delete[] zz;
}

void
CGNSMeshReader::ReadZone(Block** block, int Z)
{
    int err;

    int fn;
    err = CGNS::cg_open(mFileName.c_str(), CG_MODE_READ, &fn);
    cgns_assert(err);

    int nbases;
    err = CGNS::cg_nbases(fn, &nbases);
    cgns_assert(err);
    assert(nbases > 0);
    if (nbases > 1)
    {
        std::cerr << "Warning: there are " << nbases << " bases in the CGNS file " << mFileName
            << " but we assume the first base contains the mesh." << std::endl;
    }

    int B = 1;
    char basename[CGNS_NAME_LENGTH];
    int cell_dim, phys_dim;
    err = CGNS::cg_base_read(fn, B, basename, &cell_dim, &phys_dim);
    cgns_assert(err);

    std::cout << "cell_dim = " << cell_dim << ", phys_dim = " << phys_dim << std::endl;
    assert(cell_dim == 3);
    assert(phys_dim == 3);

    int nzones;
    err = CGNS::cg_nzones(fn, B, &nzones);
    cgns_assert(err);

    std::cout << "There are " << nzones << " zones." << std::endl;

    CGNS::CG_ZoneType_t zoneType;
    err = CGNS::cg_zone_type(fn, B, Z, &zoneType);
    cgns_assert(err);
    assert(zoneType == CGNS::CG_Structured);

    char zoneName[CGNS_NAME_LENGTH];
    CGNS::cgsize_t size[9];
    err = CGNS::cg_zone_read(fn, B, Z, zoneName, size);
    cgns_assert(err);

    std::cout << "Zone " << Z << " Name " << zoneName << std::endl;
    std::cout << "NVertices: " << size[0] << "x" << size[1] << "x" << size[2] << std::endl;
    std::cout << "NCells: " << size[3] << "x" << size[4] << "x" << size[5] << std::endl;
    int imin = 1, jmin = 1, kmin = 1, imax = size[0], jmax = size[1], kmax = size[2];

    int ncoords;
    err = CGNS::cg_ncoords(fn, B, Z, &ncoords);
    cgns_assert(err);

    std::cout << "ncoords = " << ncoords << std::endl;

    CGNS::cgsize_t rangeMin[3] = {imin, jmin, kmin};
    CGNS::cgsize_t rangeMax[3] = {imax, jmax, kmax};
    double* xx = new double[imax * jmax * kmax];
    double* yy = new double[imax * jmax * kmax];
    double* zz = new double[imax * jmax * kmax];
    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateX", CGNS::CG_RealDouble, rangeMin, rangeMax, xx);
    cgns_assert(err);
    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateY", CGNS::CG_RealDouble, rangeMin, rangeMax, yy);
    cgns_assert(err);
    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateZ", CGNS::CG_RealDouble, rangeMin, rangeMax, zz);
    cgns_assert(err);

    Close(fn);

    IndexRange meshRange(0, 0, 0, imax - 1, jmax - 1, kmax - 1);
    //*block = new Block(Z, meshRange);
    *block = Block::New(Z, meshRange);
    for (int k = 0; k < kmax; ++k)
    {
        for (int j = 0; j < jmax; ++j)
        {
            for (int i = 0; i < imax; ++i)
            {
                int ijk = k * jmax * imax + j * imax + i;
                if (i < 5 && j < 5 && k < 5)
                {
                    std::cout << "xyz: " << xx[ijk] << ' ' << yy[ijk] << ' ' << zz[ijk] << std::endl;
                }
                (*block)->XYZ()(i, j, k)[0] = xx[ijk];
                (*block)->XYZ()(i, j, k)[1] = yy[ijk];
                (*block)->XYZ()(i, j, k)[2] = zz[ijk];
            }
        }
    }
}

void
CGNSMeshReader::GetZone(int& fn, int B, int Z) const
{
    int err;

    err = CGNS::cg_open(mFileName.c_str(), CG_MODE_READ, &fn);
    cgns_assert(err);

    int nbases;
    err = CGNS::cg_nbases(fn, &nbases);
    cgns_assert(err);
    assert(nbases > 0);
    if (nbases > 1)
    {
        std::cerr << "Warning: there are " << nbases << " bases in the CGNS file " << mFileName
            << " but we assume the first base contains the mesh." << std::endl;
    }

    char basename[CGNS_NAME_LENGTH];
    int cell_dim, phys_dim;
    err = CGNS::cg_base_read(fn, B, basename, &cell_dim, &phys_dim);
    cgns_assert(err);

    std::cout << "cell_dim = " << cell_dim << ", phys_dim = " << phys_dim << std::endl;
    assert(cell_dim == 3);
    assert(phys_dim == 3);

    int nzones;
    err = CGNS::cg_nzones(fn, B, &nzones);
    cgns_assert(err);

    std::cout << "There are " << nzones << " zones." << std::endl;

    CGNS::CG_ZoneType_t zoneType;
    err = CGNS::cg_zone_type(fn, B, Z, &zoneType);
    cgns_assert(err);
    assert(zoneType == CGNS::CG_Structured);

    char zoneName[CGNS_NAME_LENGTH];
    CGNS::cgsize_t size[9];
    err = CGNS::cg_zone_read(fn, B, Z, zoneName, size);
    cgns_assert(err);

    std::cout << "Zone " << Z << " Name " << zoneName << std::endl;
    std::cout << "NVertices: " << size[0] << "x" << size[1] << "x" << size[2] << std::endl;
    std::cout << "NCells: " << size[3] << "x" << size[4] << "x" << size[5] << std::endl;
}

void
CGNSMeshReader::ReadBC(BC** bc, int Z, int BC)
{
    int err;
    int fn, B = 1;
    GetZone(fn, B, Z);

    char boconame[CGNS_NAME_LENGTH];
    CGNS::CG_BCType_t bocotype;
    CGNS::CG_PointSetType_t ptset_type;
    CGNS::cgsize_t npnts;
    int NormalIndex;
    CGNS::cgsize_t NormalListSize;
    CGNS::CG_DataType_t NormalDataType;
    int ndataset;
    err = cg_boco_info(fn, B, Z, BC, boconame, &bocotype, &ptset_type, &npnts, &NormalIndex, &NormalListSize, &NormalDataType, &ndataset);
    cgns_assert(err);
}

void
CGNSMeshReader::Close(int fn) const
{
    int err;
    err = CGNS::cg_close(fn);
    cgns_assert(err);
}

