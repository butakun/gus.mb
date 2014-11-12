// $Id: TestRamp.cpp 277 2013-06-04 01:58:51Z kato $

#include "Communicator.h"
#include "Physics.h"
#include "Roster.h"
#include "RampMesh.h"
#include "Block.h"
#include "BCSymmetry.h"
#include "BCFix.h"
#include "BCExtrapolate.h"
#include "BCCopy.h"
#include "BCViscousWall.h"
#include "ResidualEvaluator.h"
#include "Reconstructor.h"
#include "FlowModel.h"
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

int main(int argc, char** argv)
{
    int maxiter = 1000;
    //maxiter = 10;
    double cfl0 = 3.0, cflmax = 50.0, ser = 0.5;

    Communicator::Initialize(&argc, &argv);

    KOmega1988::TurbulenceModel tmodel;

    int imin = 0, jmin = 0, kmin = 0, imax = 100, jmax = 50, kmax = 1;
    IndexRange meshRange(imin, jmin, kmin, imax, jmax, kmax);
    //Block block(1, meshRange);
    //Roster::GetInstance()->RegisterBlock(0, 1);
    Block* block = Block::New(1, meshRange);

    RampMesh mesher;
    mesher.GenerateMesh(block->XYZ(), 20.0, 10.0, 1.0, 10.0, 10.0);
    block->ComputeMetrics();

    double M0 = 2.5, Rho0 = 1.2, T0 = 300.0, gamma, RGAS, C0, V0, Et;
    Physics::Initialize(1.4, Rho0, T0);
    gamma = Physics::GetInstance()->Gamma();
    RGAS = Physics::GetInstance()->RGAS();
    C0 = std::sqrt(gamma * Physics::GetInstance()->RGAS() * T0);
    V0 = M0 * C0;
    Et = RGAS / (gamma - 1.0) * T0 + 0.5 * V0 * V0;
    double rho0dim, v0dim, rhoet0dim, rho0, v0, rhoet0;
    rho0dim = Rho0;
    v0dim = V0;
    rhoet0dim = Rho0 * Et;
    rho0 = rho0dim / Physics::GetInstance()->RhoRef();
    v0 = v0dim / Physics::GetInstance()->VRef();
    rhoet0 = rhoet0dim / (Physics::GetInstance()->RhoRef() * Physics::GetInstance()->ERef());
    std::cout << "rho0dim, v0dim, rhoet0dim = " << rho0dim << ", " << v0dim << ", " << rhoet0dim << std::endl;

    double U0dim[5] = { rho0dim, v0dim, 0.0, 0.0, rhoet0dim };
    double U0[5] = { rho0, v0, 0.0, 0.0, rhoet0 };
    block->U().SetTo(U0);

    NoTurbulence turbSpec;
    tmodel.SetUT(block->UT(), block->U(), block->CellRange(), turbSpec);

    block->RegisterBC(new BCFix(IndexRange(imin, jmin, kmin, imin, jmax, kmax), I, 5, U0dim, 2, &turbSpec));
    block->RegisterBC(new BCExtrapolate(IndexRange(imax, jmin, kmin, imax, jmax, kmax), INEG));
    block->RegisterBC(new BCViscousWall(IndexRange(imin, jmin, kmin, imax, jmin, kmax), J));
    block->RegisterBC(new BCExtrapolate(IndexRange(imin, jmax, kmin, imax, jmax, kmax), JNEG));
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

    std::ofstream frms("rms.dat", std::ios::out);
    double cfl = cfl0, rms0;
    for (int iloop = 0; iloop < maxiter; ++iloop)
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
        if (iloop == 0)
        {
            rms0 = rmsTotal;
        }
        cfl = std::min(cflmax, cfl0 * std::pow(rms0 / rmsTotal, ser));
        tsEval.SetCFL(cfl);

        std::cout << iloop + 1 << '\t' << cfl << '\t' << res << std::endl;
        frms << iloop  + 1 << '\t' << rmsTotal;
        for (int l = 0; l < 5; ++l)
            frms << ' ' << rms[l];
        for (int l = 0; l < 5; ++l)
            frms << ' ' << res.MaxLocs()[l].Index;
        frms << std::endl;

        iteration.AdvanceTimeStep();
    }

    VTKWriter writer("ramp.vtk");
    writer.AddMesh(block->XYZ());
    writer.AddData(block->U(), 0, "Density", VTKWriter::SCALAR);
    writer.AddData(block->U(), 1, "Momentum", VTKWriter::VECTOR);
    writer.AddData(block->U(), 4, "TotalEnergy", VTKWriter::SCALAR);
    writer.Write();

    CGNSWriter cgnsW("ramp.cgns", false);
    int B = cgnsW.WriteBase();
    cgnsW.WriteZone(B, block->MeshRange());
    cgnsW.WriteFlowSolution(1, *block, block->U(), *Physics::GetInstance());
}

