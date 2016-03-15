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
// $Id: TestChannel.cpp 277 2013-06-04 01:58:51Z kato $

#include "Communicator.h"
#include "Physics.h"
#include "Roster.h"
#include "RectMesh.h"
#include "Block.h"
#include "BCSymmetry.h"
#include "BCFix.h"
#include "BCExtrapolate.h"
#include "BCCopy.h"
#include "BCViscousWall.h"
#include "FlowModel.h"
#include "Reconstructor.h"
#include "ResidualEvaluator.h"
#include "KOmega1988.h"
#include "SimpleExplicitIntegrator.h"
#include "LUSGSIntegrator.h"
#include "GMRESIntegrator.h"
#include "TimeStepEvaluator.h"
#include "VTKWriter.h"
#include "CGNSWriter.h"
#include "CGNSReader.h"
#include <iostream>
#include <fstream>
#include <sstream>

void Dump(const IndexRange& range, const Structured<double>& U)
{
    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                double* u = U(i, j, k);
                std::cout << IndexIJK(i, j, k) << ": " << u[0] << ", " << u[1] << ", " << u[2] << ", " << u[3] << ", " << u[4] << std::endl;
            }
        }
    }
}

void DumpVTK(const char* filename, const Block& block)
{
    VTKWriter writer(filename);
    writer.AddMesh(block.XYZ());
    writer.AddData(block.U(), 0, "Density", VTKWriter::SCALAR);
    writer.AddData(block.U(), 1, "Momentum", VTKWriter::VECTOR);
    writer.AddData(block.U(), 4, "TotalEnergy", VTKWriter::SCALAR);
    writer.Write();
}

int main(int argc, char** argv)
{
    bool unsteady = true;

    int maxiter = 1000;
    double cfl0 = 10.0, cflmax = 100.0, ser = 1.0;

    Communicator::Initialize(&argc, &argv);

    KOmega1988::TurbulenceModel tModel;

    int imin = 0, jmin = 0, kmin = 0, imax = 100, jmax = 50, kmax = 1;
    int ijetmin = imax * 0.2, ijetmax = ijetmin + 6;
    IndexRange meshRange(imin, jmin, kmin, imax, jmax, kmax);
    //Block block(1, meshRange, 2);
    Block* block = Block::New(1, meshRange, true);
    //Roster::GetInstance()->RegisterBlock(0, 1);

    RectMesh mesher;
    mesher.GenerateMesh(block->XYZ(), 1.0, 0.5, 0.01, 1.0e-3);
    block->ComputeMetrics();

    double M0 = 1.2, Rho0 = 1.2, T0 = 300.0, gamma, RGAS, C0, V0, Et;
    Physics::Initialize(1.4, Rho0, T0);
    Physics* phys = Physics::GetInstance();
    gamma = phys->Gamma();
    RGAS = phys->RGAS();
    C0 = std::sqrt(gamma * phys->RGAS() * T0);
    V0 = M0 * C0;
    Et = RGAS / (gamma - 1.0) * T0 + 0.5 * V0 * V0;
    double rho0dim, v0dim, rhoet0dim, rho0, v0, rhoet0;
    rho0dim = Rho0;
    v0dim = V0;
    rhoet0dim = Rho0 * Et;
    rho0 = rho0dim / phys->RhoRef();
    v0 = v0dim / phys->VRef();
    rhoet0 = rhoet0dim / (phys->RhoRef() * phys->ERef());
    std::cout << "rho0dim, v0dim, rhoet0dim = " << rho0dim << ", " << v0dim << ", " << rhoet0dim << std::endl;

    double U0dim[5] = { rho0dim, rho0dim * v0dim, 0.0, 0.0, rhoet0dim };
    double U0[5] = { rho0, rho0 * v0, 0.0, 0.0, rhoet0 };
    NoTurbulence turbSpec;

    double vVert = 1.0 * v0dim;
    double rhoet0Vert = rho0dim * RGAS / (gamma - 1.0) * T0 + 0.5 * rho0dim * vVert * vVert;
    double UJet[5] = { rho0dim, 0.0, rho0dim * vVert, 0.0, rhoet0Vert };

    block->U().SetTo(U0);
    //block->UT().SetTo(UT0);

    block->RegisterBC(new BCFix(IndexRange(imin, jmin, kmin, imin, jmax, kmax), I, 5, U0dim, 2, &turbSpec));
    block->RegisterBC(new BCExtrapolate(IndexRange(imax, jmin, kmin, imax, jmax, kmax), INEG));
    block->RegisterBC(new BCSymmetry(IndexRange(imin, jmin, kmin, ijetmin, jmin, kmax), J));
    BCFix* bcjet = new BCFix(IndexRange(ijetmin, jmin, kmin, ijetmax, jmin, kmax), J, 5, UJet, 2, &turbSpec);
    block->RegisterBC(bcjet);
    block->RegisterBC(new BCSymmetry(IndexRange(ijetmax, jmin, kmin, imax, jmin, kmax), J));
    block->RegisterBC(new BCSymmetry(IndexRange(imin, jmax, kmin, imax, jmax, kmax), JNEG));
    block->RegisterBC(new BCCopy(IndexRange(imin, jmin, kmin, imax, jmax, kmin), K));
    block->RegisterBC(new BCCopy(IndexRange(imin, jmin, kmax, imax, jmax, kmax), KNEG));

    block->ApplyBCs();
    block->ComputeTransportProperties();
    block->ApplyTurbBCs();
    block->ComputeTurbulentTransportProperties();

    typedef ::ResidualEvaluator<5, FirstOrderReconstructor<5> > ResEval;
    ResEval* resEval = new ResEval();
    FlowModel model;
    LUSGSIntegrator<FlowModel, ResEval> integrator(model, resEval, *block, block->U(), block->UStorage(1), block->UStorage(2));
    //SimpleExplicitIntegrator integrator(block, 0.8);
    //GMRESIntegrator<FlowModel> testIntegrator(model, block);
    int nisteps = integrator.IntegrationSteps();

    TimeStepEvaluator<FlowModel> tsEval(model, cfl0, true);
    Integrators integrators;
    integrators.push_back(&integrator);

    IterationContext iteration;
    int SUBITERATIONS = 1;
    double PERIOD, RHOVVERT0, DRHOVVERT;
    if (unsteady)
    {
        double dTime = ((1.0 / double(imax)) / v0dim) * phys->VRef() / 1.0;
        std::cout << "dTime = " << dTime << std::endl;
        iteration.SetType(IterationContext::UNSTEADY);
        iteration.SetDT(dTime);
        SUBITERATIONS = 10;
        maxiter = 100;
        cfl0 = cflmax = 100.0;
        PERIOD = dTime * 20.0;
        RHOVVERT0 = UJet[2];
        DRHOVVERT = UJet[2] * 0.5;

        CGNSReader solnReader("channel.0.cgns");
        solnReader.ReadFlowSolution(1, *block, block->U(), *phys);

        block->ApplyBCs();
        block->ComputeTransportProperties();
        block->ShiftTime();
        block->ShiftTime();
    }

    std::ofstream frms("rms.dat", std::ios::out);
    double cfl = cfl0, rms0;
    for (int iloop = 0; iloop < maxiter; ++iloop)
    {
        {
        double theta = iteration.Time() / PERIOD * 2.0 * M_PI;
        double rhovvert = RHOVVERT0 + DRHOVVERT * std::sin(theta);
        double ujet[5] = { UJet[0], UJet[1], rhovvert, UJet[3], UJet[4] };
        bcjet->SetTo(5, ujet, 2, &turbSpec);
        }

        block->ApplyBCs();
        block->ComputeTransportProperties();

        for (int isub = 0; isub < SUBITERATIONS; ++isub)
        {
            tsEval.Evaluate(integrators);
            for (int istep = 0; istep < nisteps; ++istep)
            {
                integrator.PreIntegrateStart(istep, iteration);
                integrator.PreIntegrateFinish(istep, iteration);

                integrator.Integrate(istep, iteration, block->U());
                if (istep == nisteps - 1)
                {
                    block->U().Add(1.0, integrator.DU(), block->CellRange());
                    block->ComputeTransportProperties();
                }
                integrator.PostIntegrateStart(istep, iteration);
                integrator.PostIntegrateFinish(istep, iteration);
                //Dump(IndexRange(15, 1, 0, 15, 1, 2), block->U());
            }
            Residual res = integrator.LatestResidual();
            double rms[5], rmsTotal;
            res.GetRMS(rms, rmsTotal);
            if (iloop == 0 && isub == 0)
            {
                rms0 = rmsTotal;
            }
            cfl = std::min(cflmax, cfl0 * std::pow(rms0 / rmsTotal, ser));
            tsEval.SetCFL(cfl);

            std::cout << iloop + 1 << '\t' << isub << '\t' << iteration.Time() << '\t' << cfl << '\t' << res << std::endl;
            frms << iloop  + 1 << '\t' << isub << '\t' << rmsTotal;
            for (int l = 0; l < 5; ++l)
                frms << ' ' << rms[l];
            for (int l = 0; l < 5; ++l)
                frms << ' ' << res.MaxLocs()[l].Index;
            frms << std::endl;
        }

        if (unsteady)
        {
            std::ostringstream oss;
            oss << "channel." << iloop + 1 << ".vtk";
            DumpVTK(oss.str().c_str(), *block);
        }

        block->ShiftTime();
        iteration.AdvanceTimeStep();
    }

    DumpVTK("channel.vtk", *block);

    CGNSWriter cgnsW("channel.cgns", false);
    int B = cgnsW.WriteBase();
    cgnsW.WriteZone(B, block->MeshRange());
    cgnsW.WriteFlowSolution(1, *block, block->U(), *Physics::GetInstance());
}

