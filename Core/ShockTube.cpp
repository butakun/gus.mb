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
// $Id: ShockTube.cpp 277 2013-06-04 01:58:51Z kato $

#include "Communicator.h"
#include "Physics.h"
#include "Structured.h"
#include "Block.h"
#include "RectMesh.h"
#include "ResidualEvaluator.h"
#include "BCExtrapolate.h"
#include "BCCopy.h"
#include "FlowModel.h"
#include "SimpleExplicitIntegrator.h"
#include "LUSGSIntegrator.h"
#include "TimeStepEvaluator.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

void
Dump(const Block& block, const Structured<double>& U)
{
    int j = 1, k = 1;
    for (int i = 1; i <= block.MeshRange().End.I; ++i)
    {
        std::cout << i << ": ";
        for (int l = 0; l < 5; ++l)
        {
            std::cout << U(i, j, k)[l] << " ";
        }
        std::cout << std::endl;
    }
}

void
DumpData(std::ostream& o, const Block& block, const Structured<double>& U)
{
    const Structured<double>& XYZ = block.XYZ();

    int j = 1, k = 1;
    for (int i = 1; i <= block.MeshRange().End.I; ++i)
    {
        double x = 0.5 * (XYZ(i - 1, j, k)[0] + XYZ(i, j, k)[0]);
        o << x;
        for (int l = 0; l < 5; ++l)
        {
            o << '\t' << U(i, j, k)[l];
        }
        o << std::endl;
    }
}

int
main(int argc, char** argv)
{
    Communicator::Initialize(&argc, &argv);
    Physics::Initialize();
    Physics* phys = Physics::GetInstance();

    // 10 cells in X, 1 cell in Y and Z recpectively.
    int IMAX = 100, JMAX = 1, KMAX = 1;
    //Block block(1, IndexRange(0, 0, 0, IMAX, JMAX, KMAX), 2); // 2 for two temporal storage levels
    Block* block = Block::New(1, IndexRange(0, 0, 0, IMAX, JMAX, KMAX), true);

    // 1m x 1m x 1m tube
    RectMesh mesher;
    mesher.GenerateMesh(block->XYZ(), 1.0, 1.0, 1.0);
    block->ComputeMetrics();

    // BCs
    BC* bcXi1 = new BCExtrapolate(IndexRange(0, 0, 0, 0, JMAX, KMAX), I);
    BC* bcXi2 = new BCExtrapolate(IndexRange(IMAX, 0, 0, IMAX, JMAX, KMAX), INEG);
    BC* bcEta1 = new BCCopy(IndexRange(0, 0, 0, IMAX, 0, KMAX), J);
    BC* bcEta2 = new BCCopy(IndexRange(0, JMAX, 0, IMAX, JMAX, KMAX), JNEG);
    BC* bcZeta1 = new BCCopy(IndexRange(0, 0, 0, IMAX, JMAX, 0), K);
    BC* bcZeta2 = new BCCopy(IndexRange(0, 0, KMAX, IMAX, JMAX, KMAX), KNEG);
    block->RegisterBC(bcXi1);
    block->RegisterBC(bcXi2);
    block->RegisterBC(bcEta1);
    block->RegisterBC(bcEta2);
    block->RegisterBC(bcZeta1);
    block->RegisterBC(bcZeta2);

    // Initialize the flow domain
    double T = 1.0, gamma = phys->Gamma();
    double rhoL = 1.0, pL = 1.0, rhoR = 0.1, pR = 0.1, rhoetL, rhoetR;
    rhoetL = pL / (gamma - 1.0);
    rhoetR = pR / (gamma - 1.0);
    double UL[5] = {rhoL, 0.0, 0.0, 0.0, rhoetL};
    double UR[5] = {rhoR, 0.0, 0.0, 0.0, rhoetR};

    // Acoustic speed in Region 5 (the right-most, undisturbed region)
    double c5 = std::sqrt(gamma * pR / rhoR);
    double dTime = (1.0 / double(IMAX)) / c5 * phys->VRef() / 1.0 * 0.001;
    std::cout << "dTime = " << dTime << std::endl;

    int iM = block->CellRange().End.I / 2;
    std::cout << "iM = " << iM << std::endl;

    Structured<double>& U = block->U();
    int j = block->CellRange().Start.J;
    int k = block->CellRange().Start.K;
    for (int i = block->CellRange().Start.I - 2; i <= block->CellRange().End.I + 2; ++i)
    {
        double* u;
        if (i < iM) u = UL;
        else u = UR;
        for (int l = 0; l < 5; ++l)
        {
            U(i, j, k)[l] = u[l];
        }
    }

    block->ApplyBCs();
    block->ComputeTransportProperties();

    block->ShiftTime();

    Dump(*block, block->U());

    //typedef Reconstructor<5, MinModLimiter, PrimitiveVariableCodec> Recon;
    typedef ResidualEvaluator<5, FirstOrderReconstructor<5> > ResEval;
    ResEval* resEval = new ResEval();
    FlowModel model;
    double CFL;

    int SUBITERATIONS;
#if 1
    SimpleExplicitIntegrator<FlowModel, ResEval> integrator(model, resEval, *block, block->U(), block->UStorage(1), block->UStorage(2));
    CFL = 0.8;
    TimeStepEvaluator<FlowModel> tsEval(model, CFL, false);
    IterationContext iteration;
    SUBITERATIONS = 1;
#else
    LUSGSIntegrator<FlowModel> integrator(model, *block, block->U(), block->UStorage(1), block->UStorage(2));
    CFL = 5.0;
    TimeStepEvaluator<FlowModel> tsEval(model, CFL, true);
    IterationContext iteration(IterationContext::UNSTEADY, dTime);
    SUBITERATIONS = 5;
#endif

    Integrators integrators;
    integrators.push_back(&integrator);
    int nisteps = integrator.IntegrationSteps();

    for (int i = 0; i < 50; ++i)
    {
        std::cout << "Iteration " << i << ", " << iteration << std::endl;
        for (int isub = 0; isub < SUBITERATIONS; ++isub)
        {
            if (SUBITERATIONS > 1)
                std::cout << "  Subiteration " << isub << ", " << iteration << std::endl;

            tsEval.Evaluate(integrators);

            for (int istep = 0; istep < nisteps; ++istep)
            {
                integrator.PreIntegrateStart(istep, iteration);
                integrator.PreIntegrateFinish(istep, iteration);
                integrator.Integrate(istep, iteration, U);

                if (istep == nisteps - 1)
                {
                    U.Add(1.0, integrator.DU(), block->CellRange());
                    block->ComputeTransportProperties();
                    //std::cout << "DU = " << std::endl; Dump(block, integrator.DU());
                }

                integrator.PostIntegrateStart(istep, iteration);
                integrator.PostIntegrateFinish(istep, iteration);
            }
            if (SUBITERATIONS > 1)
                std::cout << "  Residual = " << integrator.LatestResidual() << std::endl;

#if 0
                std::ostringstream oss;
                oss << std::setfill('0') << std::setw(3) << i + 1 << "_" << isub + 1 << ".dat";
                std::ofstream f(oss.str().c_str());
                DumpData(f, block, U);
#endif

            iteration.AdvanceInnerStep();
        }

        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(3) << i + 1 << ".dat";
        std::ofstream f(oss.str().c_str());
        DumpData(f, *block, U);
        Dump(*block, U);
        std::cout << oss.str() << std::endl;

        block->ShiftTime();

        iteration.AdvanceTimeStep();
    }
}

