// $Id: TestFlatPlate.cpp 277 2013-06-04 01:58:51Z kato $

#include "Communicator.h"
#include "Physics.h"
#include "Roster.h"
#include "BLMesh.h"
#include "Block.h"
#include "BCSymmetry.h"
#include "BCFix.h"
#include "BCExtrapolate.h"
#include "BCCopy.h"
#include "BCViscousWall.h"
#include "Reconstructor.h"
#include "ResidualEvaluator.h"
#include "FlowModel.h"
#include "KOmega1988.h"
#include "SimpleExplicitIntegrator.h"
#include "LUSGSIntegrator.h"
#include "TimeStepEvaluator.h"
#include "VTKWriter.h"
#include "CGNSWriter.h"
#include <iostream>
#include <fstream>
#include <iomanip>

void Dump(const IndexRange& range, const Structured<double>& U, const Structured<double>& UT, const Structured<double>& XYZ, const Structured<double>& MuK, const Structured<double>& TurMuK)
{
    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                double* u = U(i, j, k);
                double* ut = UT(i, j, k);
                double* muk = MuK(i, j, k);
                double* turmuk = TurMuK(i, j, k);
                double *p[8], xyz[3];
                p[0] = XYZ(i - 1, j - 1, k - 1);
                p[1] = XYZ(i    , j - 1, k - 1);
                p[2] = XYZ(i - 1, j    , k - 1);
                p[3] = XYZ(i    , j    , k - 1);
                p[4] = XYZ(i - 1, j - 1, k    );
                p[5] = XYZ(i    , j - 1, k    );
                p[6] = XYZ(i - 1, j    , k    );
                p[7] = XYZ(i    , j    , k    );
                xyz[0] = xyz[1] = xyz[2] = 0.0;
                for (int l = 0; l < 8; ++l)
                {
                    xyz[0] += 0.125 * p[l][0];
                    xyz[1] += 0.125 * p[l][1];
                    xyz[2] += 0.125 * p[l][2];
                }
                const size_t w = 11;
                std::cout << std::scientific << std::setprecision(4)
                << std::setw(w) << xyz[0] << std::setw(w) << xyz[1] << std::setw(w) << xyz[2]
                << std::setw(w) << u[0] << std::setw(w) << u[1] << std::setw(w) << u[2] << std::setw(w) << u[3] << std::setw(w) << u[4]
                << std::setw(w) << ut[0] << std::setw(w) << ut[1] << std::setw(w) << muk[0] << std::setw(w) << turmuk[0]
                << std::endl;
            }
        }
    }
}

void Dump(const IndexRange& r, const Structured<double>& U)
{
    for (int k = r.Start.K; k <= r.End.K; ++k)
    {
        for (int j = r.Start.J; j <= r.End.J; ++j)
        {
            for (int i = r.Start.I; i <= r.End.I; ++i)
            {
                std::cout << IndexIJK(i, j, k) << ": ";
                for (int l = 0; l < U.DOF(); ++l)
                {
                    std::cout << U(i, j, k)[l] << ", ";
                }
                std::cout << std::endl;
            }
        }
    }
}

void DumpDUDT(const IndexRange& r, const Structured<double>& DU, const Structured<double>& DT, const char* key)
{
    for (int k = r.Start.K; k <= r.End.K; ++k)
    {
        for (int j = r.Start.J; j <= r.End.J; ++j)
        {
            for (int i = r.Start.I; i <= r.End.I; ++i)
            {
                IndexIJK ii(i, j, k);
                double* du = DU(ii);
                double* dt = DT(ii);
                std::cout << key << ii << ": " << du[0] << ", " << du[1] << ", " << *dt << std::endl;
            }
        }
    }
}

void DumpProfile(std::ostream& o, const Block& block)
{
    IndexRange cr = block.CellRange();
    const Structured<double>& U = block.U();
    const Structured<double>& UT = block.UT();
    const Structured<double>& MuK = block.MuK();
    const Structured<double>& TurMuK = block.TurMuK();
    const Structured<double>& XYZ = block.XYZ();

    int i = cr.End.I - 5;
    int k = cr.Start.K;
    for (int j = cr.Start.J; j <= cr.End.J; ++j)
    {
        double y = 0.5 * (XYZ(i, j - 1, k)[1] + XYZ(i, j, k)[1]);
        double* u = U(i, j, k);
        double* ut = UT(i, j, k);
        double mu = MuK(i, j, k)[0], mut = TurMuK(i, j, k)[0];
        const size_t w = 12;
        o << std::scientific << std::setprecision(4);
        o << std::setw(w) << y << std::setw(w) << u[1] / u[0]
        << std::setw(w) << ut[0] / u[0]
        << std::setw(w) << ut[1] / u[0]
        << std::setw(w) << mu
        << std::setw(w) << mut
        << std::endl;
    }
}

int main(int argc, char** argv)
{
    bool TURBULENT = true;
    bool IMPLICIT = true;
    bool LOCALDT = true;
    int maxloop = 5000, iloopT = 10;
    double OMEGA = 1.0;

    Communicator::Initialize(&argc, &argv);

    int imin = 0, jmin = 0, kmin = 0, imax = 100, jmax = 50, kmax = 1, NLE = 10;
    IndexRange meshRange(imin, jmin, kmin, imax, jmax, kmax);
    //Block block(1, meshRange);
    //Roster::GetInstance()->RegisterBlock(0, 1);
    Block* block = Block::New(1, meshRange);

    BLMesh mesher;
    double dywall = 1.0e-5;
    //double LY = 0.005, LZ = LY / double(jmax - 1) * 0.1, LX = LZ * 20.0;
    //double LY = 0.005, LZ = dywall * 2.0, LX = 1.0;
    //double LY = 0.001, LZ = 0.001, LX = 1.0;
    double LY = 1.0, LZ = 0.0005, LX = 2.3, XLE = 0.3;
    mesher.GenerateMesh(block->XYZ(), LX, LY, LZ, NLE, XLE, dywall);
    block->ComputeMetrics();

    double M0 = 0.2, Rho0 = 0.6, T0 = 300.0, gamma, RGAS, C0, V0, Et;
    Physics::Initialize(1.4, Rho0, T0);
    gamma = Physics::GetInstance()->Gamma();
    RGAS = Physics::GetInstance()->RGAS();
    C0 = std::sqrt(gamma * Physics::GetInstance()->RGAS() * T0);
    V0 = M0 * C0;
    Et = RGAS / (gamma - 1.0) * T0 + 0.5 * V0 * V0;
    double rho0, v0, rhoet0;
    rho0 = Rho0 / Physics::GetInstance()->RhoRef();
    v0 = V0 / Physics::GetInstance()->VRef();
    rhoet0 = Rho0 * Et / (Physics::GetInstance()->RhoRef() * Physics::GetInstance()->ERef());
    std::cout << "rho0, v0, rhoet0 = " << rho0 << ", " << v0 << ", " << rhoet0 << std::endl;

    std::cout << "Reynolds number (based on the sonic speed) = " << Physics::GetInstance()->ReynoldsNumber() << std::endl;
    std::cout << "Reynolds number at the end of the plate = " << Physics::GetInstance()->ReynoldsNumber() * LX * M0 << std::endl;

    FlowModel model;
    KOmega1988::TurbulenceModel tmodel;

    double U0[5] = { rho0, rho0 * v0, 0.0, 0.0, rhoet0 };
    double U0dim[5] = { Rho0, Rho0 * V0, 0.0, 0.0, Rho0 * Et };
    block->U().SetTo(U0);

#if 0
    double Re = Physics::GetInstance()->ReynoldsNumber();
    double tke = 9.0e-9;
    double mut0 = 0.009;
    double omega = rho0 * tke / mut0;
    double rhotke = rho0 * tke;
    double rhoomega = rho0 * omega;
    double UT0[2];
#endif

    TurbulenceSpec* turbSpec;

    if (TURBULENT)
    {
        //UT0[0] = rhotke; UT0[1] = rhoomega;
        turbSpec = new TurbulenceSpecKOmegaNondimensional(9e-9, 1e-6);
    }
    else
    {
        //UT0[0] = 0.0; UT0[1] = 1.0;
        turbSpec = new NoTurbulence();
    }
    //block->UT().SetTo(UT0);
    tmodel.SetUT(block->UT(), block->U(), block->CellRange(), *turbSpec);

    block->RegisterBC(new BCFix(IndexRange(imin, jmin, kmin, imin, jmax, kmax), I, 5, U0dim, 2, turbSpec));
    block->RegisterBC(new BCExtrapolate(IndexRange(imax, jmin, kmin, imax, jmax, kmax), INEG));
#if 1
    block->RegisterBC(new BCSymmetry(IndexRange(imin, jmin, kmin, NLE, jmin, kmax), J));
    block->RegisterBC(new BCViscousWall(IndexRange(NLE, jmin, kmin, imax, jmin, kmax), J));
#else
    block->RegisterBC(new BCViscousWall(IndexRange(imin, jmin, kmin, imax, jmin, kmax), J));
#endif
    block->RegisterBC(new BCExtrapolate(IndexRange(imin, jmax, kmin, imax, jmax, kmax), JNEG));
    block->RegisterBC(new BCSymmetry(IndexRange(imin, jmin, kmin, imax, jmax, kmin), K));
    block->RegisterBC(new BCSymmetry(IndexRange(imin, jmin, kmax, imax, jmax, kmax), KNEG));

    block->ApplyBCs();
    block->ComputeTransportProperties();

    block->ApplyTurbBCs();
    block->ComputeTurbulentTransportProperties();

    typedef ::ResidualEvaluator<5, FirstOrderReconstructor<5> > ResEval;
    ResEval* resEval = new ResEval();
    KOmega1988::ResidualEvaluator* resEvalT = new KOmega1988::ResidualEvaluator();
    Integrator* integrator;
    Integrator* tintegrator;
    double cfl0, cflmax, ser, cflT0, cflTmax;
    if (IMPLICIT)
    {
        cfl0 = 5.0; cflmax = 100.0; ser = 1.0;
        cflT0 = 50.0;
        //cflT0 = 0.02;
        integrator = new LUSGSIntegrator<FlowModel, ResEval>(model, resEval, *block, block->U(), block->UStorage(1), block->UStorage(2));
        tintegrator = new LUSGSIntegrator<KOmega1988::TurbulenceModel, KOmega1988::ResidualEvaluator>(tmodel, resEvalT, *block, block->UT(), block->UT(), block->UT());
        //tintegrator = new SimpleExplicitIntegrator<KOmega1988::TurbulenceModel>(tmodel, block, block->UT());
    }
    else
    {
        cfl0 = 0.1; cflmax = 0.1; ser = 0.5;
        cflT0 = cfl0;
        integrator = new SimpleExplicitIntegrator<FlowModel, ResEval>(model, resEval, *block, block->U(), block->UStorage(1), block->UStorage(2));
        tintegrator = new SimpleExplicitIntegrator<KOmega1988::TurbulenceModel, KOmega1988::ResidualEvaluator>(tmodel, resEvalT, *block, block->UT(), block->UT(), block->UT());
    }
    int nisteps = integrator->IntegrationSteps();
    int nistepsT = tintegrator->IntegrationSteps();

    TimeStepEvaluator<FlowModel> tsEval(model, cfl0, LOCALDT);
    TimeStepEvaluator<KOmega1988::TurbulenceModel> tsEvalT(tmodel, cflT0, LOCALDT);

    Integrators integrators, tintegrators;
    integrators.push_back(integrator);
    tintegrators.push_back(tintegrator);

    IndexRange dr(imax, 0, 1, imax, jmax + 1, 1);
    IndexRange dr2(imax + 1, 0, 1, imax + 1, jmax + 1, 1);
    Dump(dr, block->U(), block->UT(), block->XYZ(), block->MuK(), block->TurMuK());

    IterationContext iteration;

    std::ofstream frms("rms.dat", std::ios::out);
    double cfl = cfl0, rms0, cflT = cflT0;
    for (int iloop = 0; iloop < maxloop; ++iloop)
    {
        if (true || iloop < iloopT)
        {
        // Flow Equations
        tsEval.Evaluate(integrators);
        for (int istep = 0; istep < nisteps; ++istep)
        {
            integrator->PreIntegrateStart(istep, iteration);
            integrator->PreIntegrateFinish(istep, iteration);

            integrator->Integrate(istep, iteration, block->U());
            if (istep == nisteps - 1)
            {
                block->U().Add(1.0, integrator->DU(), block->CellRange());
            }
            integrator->PostIntegrateStart(istep, iteration);
            integrator->PostIntegrateFinish(istep, iteration);
        }
        block->ComputeTransportProperties();
        Residual res = integrator->LatestResidual();
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
        }

        if (TURBULENT && iloop >= iloopT)
        {
            //tintegrator->DT().SetTo(integrator->DT());
            tsEvalT.Evaluate(tintegrators);

            bool copout = false;
            // Turbulence Equations
            for (int istep = 0; istep < nistepsT; ++istep)
            {
                tintegrator->PreIntegrateStart(istep, iteration);
                tintegrator->PreIntegrateFinish(istep, iteration);

                //tintegrator->SetCFL(cflT, LOCALDT);
                try
                {
                    tintegrator->Integrate(istep, iteration, block->UT());
                }
                catch (int e)
                {
                    copout = true;
                    break;
                }
                if (istep == nistepsT - 1)
                {
                    //block->UT().Add(1.0, tintegrator.DU(), block->CellRange());
                    IndexRange cr = block->CellRange();
                    for (int k = cr.Start.K; k <= cr.End.K; ++k)
                    {
                        for (int j = cr.Start.J; j <= cr.End.J; ++j)
                        {
                            for (int i = cr.Start.I; i <= cr.End.I; ++i)
                            {
                                double* ut = block->UT()(i, j, k);
                                double* dut = tintegrator->DU()(i, j, k);
                                double tmp;
                                for (int l = 0; l < 2; ++l)
                                {
                                    tmp = ut[l] + OMEGA * dut[l];
                                    ut[l] = std::max(1.0e-20, tmp);
                                }
                            }
                        }
                    }
                }
                tintegrator->PostIntegrateStart(istep, iteration);
                tintegrator->PostIntegrateFinish(istep, iteration);
            }
            if (copout)
            {
                std::cout << "Exiting" << std::endl;
                break;
            }
            block->ComputeTurbulentTransportProperties();
            //block->ApplyViscousWallBC();
            //block->FillCornerGhosts();

            Residual resT = tintegrator->LatestResidual();
            double rms[2], rmsTotal;
            resT.GetRMS(rms, rmsTotal);
            //std::cout << iloop + 1 << '\t' << "Turb RMS " << rms[0] << ' ' << rms[1] << ' ' << rmsTotal << std::endl;
            std::cout << iloop + 1 << '\t' << "Turb RMS " << resT << " " << resT.MaxLocs().size() << std::endl;
            if (false && iloop >= 0)
            {
                std::cout << "Debug = " << dr << std::endl;
                DumpDUDT(dr, tintegrator->DU(), tintegrator->DT(), "DUT, DT ");
                Dump(dr, block->U(), block->UT(), block->XYZ(), block->MuK(), block->TurMuK());
            }
        }
        if (false) {
            int dummy;
            Dump(IndexRange(2, -1, 1, 2, 30, 1), block->UT());
            std::cin >> dummy;
        }

        iteration.AdvanceTimeStep();
    }
    std::cout << "Debug = " << dr << std::endl;
    Dump(dr, block->U(), block->UT(), block->XYZ(), block->MuK(), block->TurMuK());

    std::ofstream fout("u.dat");
    DumpProfile(fout, *block);

    VTKWriter writer("flatplate.vts", VTKWriter::XML);
    writer.AddMesh(block->XYZ());
    writer.AddData(block->U(), 0, "Density", VTKWriter::SCALAR);
    writer.AddData(block->U(), 1, "Momentum", VTKWriter::VECTOR);
    writer.AddData(block->U(), 4, "TotalEnergy", VTKWriter::SCALAR);
    writer.AddData(block->UT(), 0, "RhoK", VTKWriter::SCALAR);
    writer.AddData(block->UT(), 1, "RhoOmega", VTKWriter::SCALAR);
    writer.AddData(block->MuK(), 0, "Mu", VTKWriter::SCALAR);
    writer.AddData(block->TurMuK(), 0, "TurMu", VTKWriter::SCALAR);
    writer.Write();

    //CGNSWriter cgnsWriter("flatplate.cgns");
    //cgnsWriter.WriteFlowSolution(1, block, block->U(), *Physics::GetInstance());
    //cgnsWriter.WriteTurbulenceSolution(1, block, block->U(), block->UT(), *Physics::GetInstance(), "KOmega");
}

