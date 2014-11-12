// $Id: SimpleSolver.cpp 124 2011-08-24 15:25:19Z kato $

#include "Communicator.h"
#include "Physics.h"
#include "Domain.h"
#include "DomainInfo.h"
#include "Block.h"
#include "CGNSMeshReader.h"
#include "VTKWriter.h"
#include "PLOT3DMeshWriter.h"
#include "PLOT3DSolutionWriter.h"
#include "BCSymmetry.h"
#include "BCFix.h"
#include "BCExtrapolate.h"
#include "BCViscousWall.h"
#include "BCOutletStaticPressure.h"
#include "BCInletTotal.h"
#include "Connectivity1to1.h"
#include "IndexUtils.h"
#include "FlowModel.h"
#include "SimpleExplicitIntegrator.h"
#include "LUSGSIntegrator.h"
#include "TimeStepEvaluator.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cassert>

void Dump(const Structured<double>& U, const IndexRange& range)
{
    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                std::cout << IndexIJK(i, j, k) << ": ";
                for (int l = 0; l < U.DOF(); ++l)
                {
                    std::cout << std::scientific << std::setw(12) << U(i, j, k)[l] << " ";
                }
                std::cout << std::endl;
            }
        }
    }
}

void Probe(const Block* block)
{
    return;
    const Structured<double>& U = block->U();

    if (block->ID() == 1)
    {
        std::cout << "Probing Block " << block->ID() << std::endl;
        Dump(U, IndexRange(10, 2, 2, 11, 2, 2));
    }
    else if (block->ID() == 2)
    {
        std::cout << "Probing Block " << block->ID() << std::endl;
        Dump(U, IndexRange(1, 2, 2, 2, 2, 2));
    }
    std::cout << std::endl;
}

void WritePLOT3D(const char* stem, const Blocks& blocks)
{
    for (int iblock = 0; iblock < blocks.size(); ++iblock)
    {
        Block* block = blocks[iblock];
        std::ostringstream oss;
        oss << stem << "." << std::setw(2) << std::setfill('0') << block->ID();
        std::string qname = oss.str() + ".q";
        std::string xyzname = oss.str() + ".xyz";
        PLOT3DSolutionWriter writer(qname.c_str());
        writer.Write(block->U(), block->MeshRange());

        PLOT3DMeshWriter meshWriter;
        meshWriter.Write(xyzname.c_str(), *block);
    }
}

int
main(int argc, char** argv)
{
    assert(argc > 1);

    std::ifstream f(argv[1]);
    std::string line;
    std::istringstream iss;

    Communicator::Initialize(&argc, &argv);
    Communicator* COMM = Communicator::GetInstance();
    int myrank = COMM->MyRank();

    // Mesh file
    std::getline(f, line); std::getline(f, line);
    std::string meshFileName = line;
    CGNSMeshReader meshReader(meshFileName.c_str());

    Domain domain(meshReader.Domain());

    // Scale
    std::getline(f, line); std::getline(f, line);
    iss.clear();
    iss.str(line);
    double scale = 1.0;
    iss >> scale;

    // Reference values
    std::getline(f, line); std::getline(f, line);
    iss.clear();
    iss.str(line);
    double gamma, rhoRef, TRef;
    iss >> gamma >> rhoRef >> TRef;
    Physics::Initialize(gamma, rhoRef, TRef);

    // Solver
    std::getline(f, line); std::getline(f, line);
    iss.clear();
    iss.str(line);
    int solverType = 0;
    double cfl0 = 2.0, cflmax = 1000.0, ser = 1.0;
    int local = 1;
    iss >> solverType >> cfl0 >> cflmax >> ser >> local;
    bool localTimeStepping = local != 0;

    // Iterations
    std::getline(f, line); std::getline(f, line);
    iss.clear();
    iss.str(line);
    int maxIter;
    iss >> maxIter;
    std::cout << maxIter << " iterations" << std::endl;

    // # of zones
    std::getline(f, line); std::getline(f, line);
    iss.clear();
    iss.str(line);
    int numZones;
    iss >> numZones;

    // Zones
    std::getline(f, line);
    std::map<int, int> zoneRankMap;
    for (int i = 0; i < numZones; ++i)
    {
        int zone, rank;
        std::getline(f, line);
        iss.clear(); iss.str(line);
        iss >> zone >> rank;
        zoneRankMap[zone] = rank;
    }

    // Read zones
    const DomainInfo& di = domain.GetDomainInfo();
    for (std::map<int, int>::const_iterator i = zoneRankMap.begin(); i != zoneRankMap.end(); ++i)
    {
        int Z = i->first, rank = i->second;
        const BlockInfo& bi = di.FindBlockInfo(Z);

        if (rank == myrank)
        {
            std::cout << "Reading Zone " << Z << ": Name = " << bi.Name() << ", MeshRange = " << bi.MeshRange() << std::endl;

            Block* block = new Block(bi.Zone(), bi.MeshRange());
            meshReader.ReadMesh(block->XYZ(), block->ID());
            block->XYZ() *= scale;
            block->ComputeMetrics();
            domain.RegisterLocalBlock(block);
        }
        else
        {
            domain.RegisterRemoteBlock(rank, Z);
        }
    }

    // Local zone
    for (std::map<int, int>::const_iterator i = zoneRankMap.begin(); i != zoneRankMap.end(); ++i)
    {
        int Z = i->first, rank = i->second;
        bool local = rank == myrank;

        if (!local)
        {
            std::getline(f, line);
            std::getline(f, line);
            std::getline(f, line);
            std::getline(f, line);
            int numBCs;
            iss.clear(); iss.str(line); iss >> numBCs;
            for (int ibc = 0; ibc < numBCs; ++ibc)
            {
                std::getline(f, line);
            }
            std::getline(f, line);
            std::getline(f, line);
            continue;
        }

        std::getline(f, line);
        std::getline(f, line);
        iss.clear(); iss.str(line);
        iss >> Z;
        Block* block = domain.FindLocalBlock(Z);
        const BlockInfo& bi = domain.GetDomainInfo().FindBlockInfo(Z);

        // BCs
        int numBCs;
        std::getline(f, line);
        std::getline(f, line);
        iss.clear(); iss.str(line); iss >> numBCs;

        std::cout << "Reading BCs for Zone " << Z << std::endl;
        for (int ibc = 0; ibc < numBCs; ++ibc)
        {
            std::string name, type;
            double data[5];
            std::getline(f, line);

            size_t i1, i2;
            i1 = line.find_first_of('\"', 0);
            i2 = line.find_first_of('\"', i1 + 1);
            name = line.substr(i1 + 1, i2 - i1 - 1);
            line = &line[i2 + 1];

            iss.clear(); iss.str(line);
            iss >> type;
            if (type == "Fix")
            {
                iss >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5] >> data[6];
            }
            else if (type == "ViscousWall" || type == "OutletStaticPressure")
            {
                iss >> data[0];
            }
            else if (type == "InletTotal")
            {
                iss >> data[0] >> data[1] >> data[2] >> data[3] >> data[4];
            }

            const BCInfo& bci = bi.FindBCInfoByName(name.c_str());
            Direction dir = IndexUtils::PatchDirection(bi.MeshRange(), bci.MeshRange());
            BC* bc;
            if (type == "Symmetry" || type == "InviscidWall")
            {
                bc = new BCSymmetry(bci.MeshRange(), dir);
            }
            else if (type == "Fix")
            {
                double* U0 = &data[0];
                double* UT0 = &data[5];
                bc = new BCFix(bci.MeshRange(), dir, block->U().DOF(), U0, 2, UT0);
            }
            else if (type == "Extrapolate")
            {
                bc = new BCExtrapolate(bci.MeshRange(), dir);
            }
            else if (type == "ViscousWall")
            {
                bc = new BCViscousWall(bci.MeshRange(), dir, data[0]);
            }
            else if (type == "OutletStaticPressure")
            {
                bc = new BCOutletStaticPressure(bci.MeshRange(), dir, data[0]);
            }
            else if (type == "InletTotal")
            {
                bc = new BCInletTotal(bci.MeshRange(), dir, data[0], data[1], Vector3(data[2], data[3], data[4]));
            }
            else
            {
                assert(false);
            }
            block->RegisterBC(bc);
            std::cout << "Added BC " << bc->MeshRange() << std::endl;
        }

        for (int iconn = 0; iconn < bi.GetConn1to1s().size(); ++iconn)
        {
            const Connectivity1to1Info& ci = bi.GetConn1to1s()[iconn];
            std::cout << "Conn1to1Info: " << ci << std::endl;

            Direction dir = IndexUtils::PatchDirection(bi.MeshRange(), ci.MeshRange());
            Connectivity1to1* conn = new Connectivity1to1(
                ci.MeshRange(), dir, *block, ci.Transform(), ci.DonorZone()
                );
            if (ci.GetPeriodicity() != Connectivity1to1::NONE)
            {
                conn->SetPeriodicity(ci.GetPeriodicity(), ci.GetRotationCenter(), ci.GetRotationAngle());
            }
            block->RegisterBC(conn);
            std::cout << "Added Conn1to1 " << ci.MeshRange() << std::endl;
        }

        // Initialization
        Physics* PHYS = Physics::GetInstance();
        double U0[5];
        std::getline(f, line);
        std::getline(f, line);
        iss.clear(); iss.str(line);
        iss >> U0[0] >> U0[1] >> U0[2] >> U0[3] >> U0[4];
        U0[0] /= PHYS->RhoRef();
        U0[1] /= PHYS->RhoRef() * PHYS->VRef();
        U0[2] /= PHYS->RhoRef() * PHYS->VRef();
        U0[3] /= PHYS->RhoRef() * PHYS->VRef();
        U0[4] /= PHYS->RhoRef() * PHYS->ERef();
        block->U().SetTo(U0);
        block->ApplyBCs();
        block->ComputeTransportProperties();
    }

    domain.GatherWalls();

    Blocks blocks = domain.GetLocalBlocks();

    FlowModel model;
    std::vector<Integrator*> integrators;
    for (int iblock = 0; iblock < blocks.size(); ++iblock)
    {
        Block* block = blocks[iblock];
        //integrators.push_back(new SimpleExplicitIntegrator(*block));
        integrators.push_back(new LUSGSIntegrator<FlowModel>(model, *block, block->U()));
    }

    double cfl = cfl0, rms0;
    TimeStepEvaluator<FlowModel> tsEval(model, cfl, localTimeStepping);

    COMM->Barrier();

    std::ofstream frms("rms.dat");

    for (int iloop = 0; iloop < maxIter; ++iloop)
    {
        if (COMM->MyRank() == 0)
        {
            std::cout << "Iteration " << iloop << '\t' << cfl << std::endl;
        }

        tsEval.Evaluate(integrators);

        int nisteps = integrators[0]->IntegrationSteps();
        for (int ii = 0; ii < nisteps; ++ii)
        {
            // Pre-Integrate Phase
            for (int iblock = 0; iblock < blocks.size(); ++iblock)
            {
                integrators[iblock]->PreIntegrateStart(ii);
            }
            for (int iblock = 0; iblock < blocks.size(); ++iblock)
            {
                integrators[iblock]->PreIntegrateFinish(ii);
            }

            // Integrate Phase
            for (int iblock = 0; iblock < blocks.size(); ++iblock)
            {
                Integrator* integrator = integrators[iblock];
                integrator->Integrate(ii, blocks[iblock]->U());
            }
            if (ii == nisteps - 1)
            {
                for (int iblock = 0; iblock < blocks.size(); ++iblock)
                {
                    Block* block = blocks[iblock];
                    Integrator* integrator = integrators[iblock];
                    block->U().Add(1.0, integrator->DU(), block->CellRange());
                    block->ComputeTransportProperties();
                }
            }

            // Post-Integrate Phase
            for (int iblock = 0; iblock < blocks.size(); ++iblock)
            {
                integrators[iblock]->PostIntegrateStart(ii);
            }
            for (int iblock = 0; iblock < blocks.size(); ++iblock)
            {
                integrators[iblock]->PostIntegrateFinish(ii);
                Probe(blocks[iblock]);
            }
            COMM->Barrier();
        }

        // Gather and accumulate the residuals from all the local blocks
        Residual res = integrators[0]->LatestResidual();
        for (int iblock = 1; iblock < blocks.size(); ++iblock)
        {
            res += integrators[iblock]->LatestResidual();
        }
        // Then gather and accumulate the residuals from all the blocks globally
        res.Allgather();

        // Dump the residual 
        double rms[5], rmsTotal;
        res.GetRMS(rms, rmsTotal);
        COMM->Console() << iloop;
        for (int l = 0; l < 5; ++l)
        {
            COMM->Console() << '\t' << rms[l];
        }
        COMM->Console() << '\t' << rmsTotal;
        for (int l = 0; l < 5; ++l)
        {
            COMM->Console() << '\t' << "(" << res.MaxLocs()[l].BlockID << ", " << res.MaxLocs()[l].Index << ")";
        }
        COMM->Console() << std::endl;

        if (COMM->MyRank() == 0)
        {
            for (int l = 0; l < 5; ++l)
            {
                std::cout << '\t' << rms[l];
            }
            std::cout << '\t' << rmsTotal << std::endl;
            frms << iloop << '\t' << rmsTotal << std::endl;
        }

        if (iloop == 0)
        {
            rms0 = rmsTotal;
        }
        cfl = std::max(cfl0, std::min(cflmax, cfl0 * std::pow(rms0 / rmsTotal, ser)));
        tsEval.SetCFL(cfl);

        if (iloop % 1000 == 0)
        {
            std::ostringstream oss;
            oss << "out." << std::setw(6) << std::setfill('0') << iloop;
            WritePLOT3D(oss.str().c_str(), blocks);
        }
    }
    std::cout << "DONE" << std::endl;

    if (true || COMM->MyRank() == 0)
    {
    for (int iblock = 0; iblock < blocks.size(); ++iblock)
    {
        Block* block = blocks[iblock];
        std::ostringstream oss;
        oss << "out." << std::setw(2) << std::setfill('0') << block->ID() << ".vtk";
        VTKWriter writer(oss.str().c_str());
        writer.AddMesh(block->XYZ());
        writer.AddData(block->U(), 0, "Density", VTKWriter::SCALAR);
        writer.AddData(block->U(), 4, "TotalEnergy", VTKWriter::SCALAR);
        writer.AddData(block->U(), 1, "Momentum", VTKWriter::VECTOR);
        writer.Write();
    }
    }

#if 0
    for (int iblock = 0; iblock < blocks.size(); ++iblock)
    {
        Block* block = blocks[iblock];
        std::ostringstream oss;
        oss << "out." << std::setw(2) << std::setfill('0') << block->ID();
        std::string qname = oss.str() + ".q";
        std::string xyzname = oss.str() + ".xyz";
        PLOT3DSolutionWriter writer(qname.c_str());
        writer.Write(block->U(), block->MeshRange());

        PLOT3DMeshWriter meshWriter;
        meshWriter.Write(xyzname.c_str(), *block);
    }
#else
    WritePLOT3D("final", blocks);
#endif

    COMM->Finalize();

    return EXIT_SUCCESS;
}

