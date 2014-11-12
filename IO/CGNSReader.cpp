// $Id: CGNSReader.cpp 321 2014-08-26 09:52:44Z kato $

#include "CGNSReader.h"
#include "Block.h"
#include <cassert>
namespace CGNS {
#include <cgnslib.h>
}

const size_t CGNS_NAME_LENGTH = 33;

void cgns_assert(int err)
{
    if (err == 0)
        return;
    std::cerr << CGNS::cg_get_error() << std::endl;
    assert(false);
}

CGNSReader::CGNSReader(const char* filename)
:   fn(-1)
{
    int err;

    err = CGNS::cg_open(filename, CG_MODE_READ, &this->fn);
    cgns_assert(err);
}

CGNSReader::~CGNSReader()
{
    int err;
    err = CGNS::cg_close(fn);
    cgns_assert(err);
}

CGNSStructure
CGNSReader::ReadStructure() const
{
    CGNSStructure s;
    int err;

    int nbases;
    err = CGNS::cg_nbases(fn, &nbases);
    cgns_assert(err);

    for (int B = 1; B <= nbases; ++B)
    {
        char basename[CGNS_NAME_LENGTH];
        int cell_dim, phys_dim;
        err = CGNS::cg_base_read(fn, B, basename, &cell_dim, &phys_dim);
        cgns_assert(err);

        CGNSStructure::Base base(B, basename, cell_dim, phys_dim);

        int nzones;
        err = CGNS::cg_nzones(fn, B, &nzones);
        cgns_assert(err);

        for (int Z = 1; Z <= nzones; ++Z)
        {
            char zonename[CGNS_NAME_LENGTH];
            CGNS::cgsize_t size[9];
            err = CGNS::cg_zone_read(fn, B, Z, zonename, size);
            cgns_assert(err);

            CGNSStructure::Zone zone(Z, zonename, IndexIJK(size[0], size[1], size[2]), IndexIJK(size[3], size[4], size[5]));
            base.Zones.push_back(zone);
        }

        s.Bases.push_back(base);
    }

    return s;
}

void
CGNSReader::ReadMesh(int zone, Structured<double>& XYZ, std::string& zoneName, const IndexRange rangeToRead)
{
    int B = 1;
    int err;

    char zoneName_[CGNS_NAME_LENGTH];
    CGNS::cgsize_t size[9];
    err = CGNS::cg_zone_read(this->fn, B, zone, zoneName_, size);
    zoneName = zoneName_;

    IndexRange mr = rangeToRead;
    if (rangeToRead == IndexRange(0, 0, 0, 0, 0, 0))
    {
        mr = IndexRange(0, 0, 0, size[0] -1, size[1] - 1, size[2] - 1);
    }

    CGNS::cgsize_t rangeMin[3] = {mr.Start.I + 1, mr.Start.J + 1, mr.Start.K + 1};
    CGNS::cgsize_t rangeMax[3] = {mr.End.I + 1, mr.End.J + 1, mr.End.K + 1};

    // Temporary storage for the coord values
    double* xx = new double[mr.Count()];
    double* yy = new double[mr.Count()];
    double* zz = new double[mr.Count()];

    err = CGNS::cg_coord_read(fn, B, zone, "CoordinateX", CGNS::CG_RealDouble, rangeMin, rangeMax, xx);
    cgns_assert(err);
    err = CGNS::cg_coord_read(fn, B, zone, "CoordinateY", CGNS::CG_RealDouble, rangeMin, rangeMax, yy);
    cgns_assert(err);
    err = CGNS::cg_coord_read(fn, B, zone, "CoordinateZ", CGNS::CG_RealDouble, rangeMin, rangeMax, zz);
    cgns_assert(err);

    Structured<double> XX(xx, 1, mr);
    Structured<double> YY(yy, 1, mr);
    Structured<double> ZZ(zz, 1, mr);

    XYZ.Allocate(3, mr);
    for (int k = mr.Start.K; k <= mr.End.K; ++k)
    {
        for (int j = mr.Start.J; j <= mr.End.J; ++j)
        {
            for (int i = mr.Start.I; i <= mr.End.I; ++i)
            {
                double* xyz = XYZ(i, j, k);
                xyz[0] = *XX(i, j, k);
                xyz[1] = *YY(i, j, k);
                xyz[2] = *ZZ(i, j, k);
            }
        }
    }

    delete[] xx;
    delete[] yy;
    delete[] zz;
}

void
CGNSReader::ReadFlowSolution(int zone, const Block& block, Structured<double>& U, const Physics& phys)
{
    int err;
    int B = 1;

    char solname[CGNS_NAME_LENGTH];
    CGNS::CG_GridLocation_t location;
    int S = 1;
    err = CGNS::cg_sol_info(this->fn, B, zone, S, solname, &location);
    cgns_assert(err);

    std::string solName(solname);
    assert(solName == "FlowSolution");

    IndexRange cr = block.CellRange();

    double* rho = new double[cr.Count()];
    double* rhou = new double[cr.Count()];
    double* rhov = new double[cr.Count()];
    double* rhow = new double[cr.Count()];
    double* rhoet = new double[cr.Count()];
    Structured<double> Rho(rho, 1, cr);
    Structured<double> RhoU(rhou, 1, cr);
    Structured<double> RhoV(rhov, 1, cr);
    Structured<double> RhoW(rhow, 1, cr);
    Structured<double> RhoEt(rhoet, 1, cr);

    CGNS::cgsize_t rangeMin[3], rangeMax[3];
    // rangeMin/Max are 1-bias (CGNS convention), cr is 0-bias but our first internal cell is cr.Start + 1, so...
    rangeMin[0] = cr.Start.I; rangeMin[1] = cr.Start.J; rangeMin[2] = cr.Start.K;
    rangeMax[0] = cr.End.I; rangeMax[1] = cr.End.J; rangeMax[2] = cr.End.K;

    err = CGNS::cg_field_read(this->fn, B, zone, S, "Density", CGNS::CG_RealDouble, rangeMin, rangeMax, rho);
    cgns_assert(err);
    err = CGNS::cg_field_read(this->fn, B, zone, S, "MomentumX", CGNS::CG_RealDouble, rangeMin, rangeMax, rhou);
    cgns_assert(err);
    err = CGNS::cg_field_read(this->fn, B, zone, S, "MomentumY", CGNS::CG_RealDouble, rangeMin, rangeMax, rhov);
    cgns_assert(err);
    err = CGNS::cg_field_read(this->fn, B, zone, S, "MomentumZ", CGNS::CG_RealDouble, rangeMin, rangeMax, rhow);
    cgns_assert(err);
    err = CGNS::cg_field_read(this->fn, B, zone, S, "EnergyStagnationDensity", CGNS::CG_RealDouble, rangeMin, rangeMax, rhoet);
    cgns_assert(err);

    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* u = U(i, j, k);
                u[0] = *Rho(i, j, k) / phys.RhoRef();
                u[1] = *RhoU(i, j, k) / (phys.RhoRef() * phys.VRef());
                u[2] = *RhoV(i, j, k) / (phys.RhoRef() * phys.VRef());
                u[3] = *RhoW(i, j, k) / (phys.RhoRef() * phys.VRef());
                u[4] = *RhoEt(i, j, k) / (phys.RhoRef() * phys.ERef());
            }
        }
    }

    delete[] rho;
    delete[] rhou;
    delete[] rhov;
    delete[] rhow;
    delete[] rhoet;
}

void
CGNSReader::ReadTurbulenceSolution(int zone, const Block& block, const Structured<double>& U, Structured<double>& UT, const Physics& phys, const std::string& model)
{
    int err;
    int B = 1;

    char solname[CGNS_NAME_LENGTH];
    CGNS::CG_GridLocation_t location;
    int S = 1;
    err = CGNS::cg_sol_info(this->fn, B, zone, S, solname, &location);
    cgns_assert(err);

    std::string solName(solname);
    assert(solName == "FlowSolution");

    assert(model == "KOmega");

    IndexRange cr = block.CellRange();

    double* tke = new double[cr.Count()];
    double* omega = new double[cr.Count()];
    Structured<double> TKE(tke, 1, cr);
    Structured<double> Omega(omega, 1, cr);

    CGNS::cgsize_t rangeMin[3], rangeMax[3];
    // rangeMin/Max are 1-bias (CGNS convention), cr is 0-bias but our first internal cell is cr.Start + 1, so...
    rangeMin[0] = cr.Start.I; rangeMin[1] = cr.Start.J; rangeMin[2] = cr.Start.K;
    rangeMax[0] = cr.End.I; rangeMax[1] = cr.End.J; rangeMax[2] = cr.End.K;

    err = CGNS::cg_field_read(this->fn, B, zone, S, "TurbulentEnergyKinetic", CGNS::CG_RealDouble, rangeMin, rangeMax, tke);
    cgns_assert(err);
    err = CGNS::cg_field_read(this->fn, B, zone, S, "TurbulentDissipationRate", CGNS::CG_RealDouble, rangeMin, rangeMax, omega);
    cgns_assert(err);

    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* u = U(i, j, k);
                double* ut = UT(i, j, k);
                ut[0] = *TKE(i, j, k) / phys.TKERef() * u[0];
                ut[1] = *Omega(i, j, k) / phys.OmegaRef() * u[0];
            }
        }
    }

    delete[] tke;
    delete[] omega;
}

