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

#include "CGNSWriter.h"
#include "CGNSReader.h"
#include "Structured.h"
#include <sstream>
#include <cassert>
namespace CGNS {
#include <cgnslib.h>
}

const size_t CGNS_NAME_LENGTH = 33;

static void cgns_assert(int err)
{
    if (err == 0)
        return;
    std::cerr << CGNS::cg_get_error() << std::endl;
    assert(false);
}

CGNSWriter::CGNSWriter(const char* filename, bool modify)
{
    int err;

    int mode = modify ? CG_MODE_MODIFY : CG_MODE_WRITE;
    err = CGNS::cg_open(filename, mode, &this->fn);
    cgns_assert(err);
}

CGNSWriter::~CGNSWriter()
{
    int err;
    err = CGNS::cg_close(fn);
    cgns_assert(err);
}

void
CGNSWriter::WriteStructure(const CGNSStructure& s, const Physics& phys)
{
    int err;

    for (CGNSStructure::Bases_t::const_iterator ib = s.Bases.begin(); ib != s.Bases.end(); ++ib)
    {
        const CGNSStructure::Base& base = *ib;

        int B;
        err = CGNS::cg_base_write(fn, base.Name.c_str(), base.CellDim, base.PhysDim, &B);
        cgns_assert(err);
        if (B != base.B)
            std::cerr << "CGNSWriter: base number B has changed." << std::endl;

        for (CGNSStructure::Zones_t::const_iterator iz = base.Zones.begin(); iz != base.Zones.end(); ++iz)
        {
            const CGNSStructure::Zone& zone = *iz;

            int Z;
            CGNS::cgsize_t size[9];
            size[0] = zone.VertexSize.I;
            size[1] = zone.VertexSize.J;
            size[2] = zone.VertexSize.K;
            size[3] = zone.CellSize.I;
            size[4] = zone.CellSize.J;
            size[5] = zone.CellSize.K;
            size[6] = size[7] = size[8] = 0;
            err = CGNS::cg_zone_write(fn, B, zone.Name.c_str(), size, CGNS::CG_Structured, &Z);
            cgns_assert(err);
            if (Z != zone.Z)
                std::cerr << "CGNSWriter: zone number Z has changed." << std::endl;

            RigidBodyMotion* rbm = zone.Motion;
            if (rbm != NULL)
            {
                RotationalMotion* rm = dynamic_cast<RotationalMotion*>(rbm);
                if (rm != NULL)
                {
                    int R;
                    err = CGNS::cg_rigid_motion_write(fn, B, Z, "RigidMotion", CGNS::CG_ConstantRate, &R);
                    cgns_assert(err);

                    std::ostringstream oss;
                    oss << "/" << base.Name << "/" << zone.Name << "/" << "RigidMotion";
                    err = CGNS::cg_gopath(fn, oss.str().c_str());
                    cgns_assert(err);

                    Vector3 o = rm->Origin() * phys.LRef();
                    Vector3 omega = rm->AngularVelocity() / phys.TimeRef();
                    double origin[6];
                    origin[0] = o.X(); origin[1] = o.Y(); origin[2] = o.Z();
                    origin[3] = o.X(); origin[4] = o.Y(); origin[5] = o.Z();
                    CGNS::cgsize_t dims[2] = { 3, 2 };
                    err = CGNS::cg_array_write("OriginLocation", CGNS::CG_RealDouble, 2, dims, origin);
                    cgns_assert(err);

                    double rotation[3];
                    rotation[0] = omega.X(); rotation[1] = omega.Y(); rotation[2] = omega.Z();
                    dims[0] = 3;
                    err = CGNS::cg_array_write("RigidRotationRate", CGNS::CG_RealDouble, 1, dims, rotation);
                    cgns_assert(err);
                }
            }
        }
    }
}

int
CGNSWriter::WriteBase()
{
    int B;
    int err;
    err = CGNS::cg_base_write(this->fn, "base", 3, 3, &B);
    cgns_assert(err);
    return B;
}

int
CGNSWriter::WriteZone(int B, const IndexRange& mr, const char* name)
{
    std::string zoneName = "FakeName";
    if (name != NULL)
    {
        zoneName = name;
    }

    CGNS::cgsize_t size[9];
    size[0] = mr.Size().I; size[1] = mr.Size().J; size[2] = mr.Size().K;
    size[3] = mr.Size().I - 1; size[4] = mr.Size().J - 1; size[5] = mr.Size().K - 1;
    size[6] = size[7] = size[8] = 0;
    int Z;
    int err;
    err = CGNS::cg_zone_write(this->fn, B, zoneName.c_str(), size, CGNS::CG_Structured, &Z);
    cgns_assert(err);

    return Z;
}

void
CGNSWriter::WriteFlowSolution(int zone, const Block& block, const Structured<double>& U, const Physics& phys)
{
    int err;
    int B = 1;
    int S;

    IndexRange mr = block.MeshRange();
    IndexRange cr = block.CellRange();

    err = CGNS::cg_sol_write(this->fn, B, zone, "FlowSolution", CGNS::CG_CellCenter, &S);
    cgns_assert(err);

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

    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* u = U(i, j, k);
                *Rho(i, j, k) = u[0] * phys.RhoRef();
                *RhoU(i, j, k) = u[1] * phys.RhoRef() * phys.VRef();
                *RhoV(i, j, k) = u[2] * phys.RhoRef() * phys.VRef();
                *RhoW(i, j, k) = u[3] * phys.RhoRef() * phys.VRef();
                *RhoEt(i, j, k) = u[4] * phys.RhoRef() * phys.ERef();
            }
        }
    }

    // FIXME: MomentumXYZ must be in inertial frame of reference, or use RotationalMomentumXYZ
    int f1, f2, f3, f4, f5;
    err = CGNS::cg_field_write(fn, B, zone, S, CGNS::CG_RealDouble, "Density", rho, &f1);
    cgns_assert(err);
    err = CGNS::cg_field_write(fn, B, zone, S, CGNS::CG_RealDouble, "MomentumX", rhou, &f2);
    cgns_assert(err);
    err = CGNS::cg_field_write(fn, B, zone, S, CGNS::CG_RealDouble, "MomentumY", rhov, &f3);
    cgns_assert(err);
    err = CGNS::cg_field_write(fn, B, zone, S, CGNS::CG_RealDouble, "MomentumZ", rhow, &f4);
    cgns_assert(err);
    err = CGNS::cg_field_write(fn, B, zone, S, CGNS::CG_RealDouble, "EnergyStagnationDensity", rhoet, &f5);
    cgns_assert(err);

    std::cout << "CGNSWriter: f1, f2, f3, f4, f5 = " << f1 << ", " << f2 << ", " << f3 << ", " << f4 << ", " << f5 << std::endl;

    delete[] rho;
    delete[] rhou;
    delete[] rhov;
    delete[] rhow;
    delete[] rhoet;
}

void
CGNSWriter::WriteTurbulenceSolution(int zone, const Block& block, const Structured<double>& U, const Structured<double>& UT, const Physics& phys, const std::string& model)
{
    // Assumes the primary flow solution (see WriteFlowSolution) has been written.

    int err;
    int B = 1;
    int S = 1;
    char solname[CGNS_NAME_LENGTH];
    CGNS::CG_GridLocation_t location;

    IndexRange mr = block.MeshRange();
    IndexRange cr = block.CellRange();

    err = CGNS::cg_sol_info(this->fn, B, zone, S, solname, &location);
    cgns_assert(err);

    assert(std::string(solname) == "FlowSolution");

    assert(model == "KOmega");

    double* turbk = new double[cr.Count()];
    double* turbomega = new double[cr.Count()];
    Structured<double> TurbK(turbk, 1, cr);
    Structured<double> TurbOmega(turbomega, 1, cr);

    double CC = phys.VRef() * phys.VRef();
    double OmegaRef = phys.RhoRef() * CC / phys.MuRef();

    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* u = U(i, j, k);
                double* ut = UT(i, j, k);
                *TurbK(i, j, k) = ut[0] / u[0] * phys.TKERef();
                *TurbOmega(i, j, k) = ut[1] / u[0] * phys.OmegaRef();
            }
        }
    }

    // FIXME: MomentumXYZ must be in inertial frame of reference
    int f1, f2;
    err = CGNS::cg_field_write(fn, B, zone, S, CGNS::CG_RealDouble, "TurbulentEnergyKinetic", turbk, &f1);
    cgns_assert(err);
    err = CGNS::cg_field_write(fn, B, zone, S, CGNS::CG_RealDouble, "TurbulentDissipationRate", turbomega, &f2);
    cgns_assert(err);

    delete[] turbk;
    delete[] turbomega;
}

