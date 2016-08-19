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

#include "Communicator.h"
#include "Profiler.h"
#include "Roster.h"
#include "Physics.h"
#include "Block.h"
#include "VTKWriter.h"
#include "PLOT3DMeshWriter.h"
#include "PLOT3DSolutionWriter.h"
#include "BCSymmetry.h"
#include "BCFix.h"
#include "BCExtrapolate.h"
#include "BCViscousWall.h"
#include "BCRotatingViscousWall.h"
#include "BCOutletStaticPressure.h"
#include "BCInletTotal.h"
#include "BCEnforce2D.h"
#include "BCForcePull.h"
#include "Connectivity1to1.h"
#include "ConnectivityAbutting.h"
#include "IndexUtils.h"
#include "FlowModel.h"
#include "ResidualEvaluator.h"
#include "Reconstructor.h"
#include "KOmega1988.h"
#include "IterationContext.h"
#include "SimpleExplicitIntegrator.h"
#include "LUSGSIntegrator.h"
#include "TimeStepEvaluator.h"
#include "StructuredDataExchanger.h"
#include "SimplePlanarAbuttingInterface.h"
#include "InterfaceDataExchanger.h"
#include "CGNSReader.h"
#include "CGNSWriter.h"
#include "FieldInitializer.h"
#include "Matrix33.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>

class BlockOffset
{
public:
    class Rotate
    {
    public:
        Matrix33 R; 
        Rotate(const Vector3& v) { R = Matrix33::RotationMatrix(v); }
        void operator () (double* xyz) const
        {
            R.Multiply(xyz);
        }
    };

    enum OffsetType { NONE, ROTATION, TRANSLATION } Type;
    Vector3 Offset;

    BlockOffset(OffsetType type, const Vector3& offset) : Type(type), Offset(offset) {}
    BlockOffset() : Type(NONE), Offset(0.0, 0.0, 0.0) {}

    void ApplyOffset(Structured<double>& XYZ) const
    {
        if (Type == NONE)
        {
            return;
        }
        else if (Type == ROTATION)
        {
            std::ostream& LOG = Communicator::GetInstance()->Console();
            Rotate rotate(Offset);
            LOG << " BlockOffset: rotation matrix = " << std::endl << rotate.R << std::endl;
            XYZ.Apply(rotate);
        }
        else
        {
            assert(false);
        }
    }
};

class BlockDef
{
public:
    BlockDef() {}
    BlockDef(int b, int rank, int z, const IndexRange& meshRange, const BlockOffset& offset)
    : BlockID(b), Rank(rank), Z(z), MeshRange(meshRange), Offset(offset), LocalBlock(NULL) {}

    int BlockID;
    int Rank;
    int Z;
    IndexRange MeshRange;
    BlockOffset Offset;
    Block* LocalBlock;
};

typedef std::map<int, BlockDef> BlockRankMap;

class BlockRange
{
public:
    int BlockID;
    IndexRange Range;
    BlockRange(int blockID, const IndexRange& range)
    : BlockID(blockID), Range(range)
    {}
};

std::ostream& operator << (std::ostream& o, const BlockRange& br)
{
    o << br.BlockID << ":" << br.Range;
    return o;
}

typedef std::vector<BlockRange> BlockRanges;

void splash(std::ostream& o)
{
    o << "gus.mb  Copyright (C) 2016  Hiromasa Kato <hiromasa at gmail.com>" << std::endl
      << "This program comes with ABSOLUTELY NO WARRANTY." << std::endl
      << "This is free software, and you are welcome to redistribute it" << std::endl
      << "under certain conditions." << std::endl;
}

void
WriteSolution(const char* filename, Blocks& blocks, const BlockRankMap& blockRankMap)
{
    Communicator* comm = Communicator::GetInstance();

    if (comm->MyRank() == 0)
    {
        CGNSWriter writer(filename, false);
        int B;
        B = writer.WriteBase();
        for (int ib = 1; ib <= blockRankMap.size(); ++ib)
        {
            const VirtualBlock* block = Roster::GetInstance()->GetBlock(ib);
            std::ostringstream oss;
            oss << block->ID();
            writer.WriteZone(B, block->MeshRange(), oss.str().c_str());
        }
    }

    for (int rank = 0; rank < comm->Size(); ++rank)
    {
        if (rank == comm->MyRank())
        {
            CGNSWriter writer(filename);
            for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
            {
                const Block* block = blocks[iblock];
                BlockRankMap::const_iterator itor = blockRankMap.find(block->ID());
                const BlockDef& bdef = itor->second;
                std::cout << "Writing Block " << block->ID() << " solution" << std::endl;
                writer.WriteFlowSolution(bdef.Z, *block, block->U(), *Physics::GetInstance());
            }
        }

        comm->Barrier();
    }
}

void Enforce2D(Structured<double>& U)
{
    for (IndexIterator it(U.GetRange()); !it.IsEnd(); it.Advance())
    {
        IndexIJK ijk = it.Index();
        double* u = U(ijk);
        double P = std::sqrt(u[1] * u[1] + u[2] * u[2] + u[3] * u[3]); // Momentum
        double e, e1, e2;
        e = std::sqrt(u[1] * u[1] + u[2] * u[2]);
        e1 = u[1] / e;
        e2 = u[2] / e;
        u[1] = P * e1;
        u[2] = P * e2;
        u[3] = 0.0;
    }
}

void Enforce2D(Blocks& blocks)
{
    for (Blocks::iterator i = blocks.begin(); i != blocks.end(); ++i)
    {
        Enforce2D((*i)->U());
    }
}

void Initialize(const FlowModel& model, const Block& block, Structured<double>& U, double* UGlobal)
{
    int dof = model.DOF();
    double* ULocal = new double[dof];

    for (IndexIterator itor(block.CellRange()); !itor.IsEnd(); itor.Advance())
    {
        IndexIJK ijk = itor.Index();
        model.FromGlobalToLocal(ULocal, UGlobal, block, ijk);
        double* UU = U(ijk);
        for (int l = 0; l < dof; ++l)
        {
            UU[l] = ULocal[l];
        }
    }

    delete[] ULocal;
}

void Probe(std::ostream& out, const Block& block, const IndexRange& range)
{
    const Structured<double>& XYZ = block.XYZ();
    const Structured<double>& U = block.U();
    const Structured<double>& UT = block.UT();
    const Structured<double>& MuK = block.MuK();
    const Structured<double>& TurMuK = block.TurMuK();
    const Structured<double>& Rad = block.Radius();

    out << "PROBE_BEGIN: U[0,1,2,3,4], UT[0,1], TurMuK[0], MuK[0], rad[0,1,2,3]" << std::endl;
    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                double* rad = Rad(i, j, k);
                Vector3 pC(rad[0], rad[1], rad[2]);

                out << block.ID() << ":" << IndexIJK(i, j, k) << ": ";
                out << std::scientific << std::setw(12) << std::setprecision(4);
                out << U(i, j, k)[0] << " " << U(i, j, k)[1] << " " << U(i, j, k)[2] << " " << U(i, j, k)[3] << " " << U(i, j, k)[4] << " ";
                out << UT(i, j, k)[0] << " " << UT(i, j, k)[1] << " " << TurMuK(i, j, k)[0] << " " << MuK(i, j, k)[0] << " ";
                out << rad[0] << " " << rad[1] << " " << rad[2] << " " << rad[3];
                out << std::endl;
            }
        }
    }
    out << "PROBE_END" << std::endl;
}

void Probe(std::ostream& out, const Blocks& blocks, const BlockRanges& ranges)
{
    for (BlockRanges::const_iterator i = ranges.begin(); i != ranges.end(); ++i)
    {
        const BlockRange& br = *i;
        const Block* block = NULL;
        for (Blocks::const_iterator ii = blocks.begin(); ii != blocks.end(); ++ii)
        {
            if ((*ii)->ID() == br.BlockID)
            {
                block = *ii;
                break;
            }
        }
        if (block != NULL)
        {
            Probe(out, *block, br.Range);
        }
    }
}

void DumpData(std::ostream& out, const Structured<double>& data, const IndexRange& range)
{
    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                out << "DUMP:";
                out << IndexIJK(i, j, k) << ": ";
                out << std::scientific << std::setw(12) << std::setprecision(4);
                for (int l = 0; l < data.DOF(); ++l)
                {
                    out << data(i, j, k)[l] << " ";
                }
                out << std::endl;
            }
        }
    }
}

void ReadBlockRanges(BlockRanges& ranges, const char* filename)
{
    std::ifstream f(filename);
    std::string line;

    while (std::getline(f, line))
    {
        std::istringstream  iss(line);
        char c;
        iss >> c;
        if (c == '#')
            continue;

        int blockID, i1, j1, k1, i2, j2, k2;
        iss.clear();
        iss.str(line);
        iss >> blockID >> i1 >> j1 >> k1 >> i2 >> j2 >> k2;
        if (iss)
        {
            BlockRange br(blockID, IndexRange(i1, j1, k1, i2, j2, k2));
            ranges.push_back(br);
            std::ostream& CONSOLE = Communicator::GetInstance()->Console();
            CONSOLE << "Debug output enabled for " << br << std::endl;
        }
    }

    f.close();
}

TurbulenceSpec* ReadTurbulenceSpec(std::istringstream& iss)
{
    std::string type;
    iss >> type;
    if (type == "IntensityViscosityRatio")
    {
        double I, viscRatio, rho, U, T;
        iss >> I >> viscRatio >> rho >> U >> T;
        return new TurbulenceSpecIntensityViscosityRatio(I, viscRatio, rho, U, T);
    }
    else if (type == "NoTurbulence")
    {
        return new NoTurbulence();
    }
    else if (type == "KOmegaNondimensional")
    {
        double tke, omega;
        iss >> tke >> omega;
        return new TurbulenceSpecKOmegaNondimensional(tke, omega);
    }

    return NULL;
}

std::string UnQuote(const std::string& quoted)
{
    size_t i1 = quoted.find_first_of('\"', 0);
    size_t i2 = quoted.find_first_of('\"', i1 + 1);
    std::string unquoted = quoted.substr(i1 + 1, i2 - i1 - 1);
    return unquoted;
}

void WriteSolutionCGNS(const char* solutionStem, const Blocks& blocks, const BlockRankMap& blockRankMap, const CGNSStructure& cgnsStruct)
{
    Communicator* COMM = Communicator::GetInstance();
    Profiler* PROF = Profiler::GetInstance();

    PROF->CheckPoint();

    // Create a new CGNS file and initialize the topology.
    std::string solutionFileName = solutionStem;
    solutionFileName += ".cgns";
    if (COMM->MyRank() == 0)
    {
        CGNSWriter* writer = new CGNSWriter(solutionFileName.c_str(), false);
        writer->WriteStructure(cgnsStruct);
        delete writer; // deleting means closing the file.
    }
    COMM->Barrier();

    for (int rank = 0; rank < COMM->Size(); ++rank)
    {
        if (rank == COMM->MyRank())
        {
            CGNSWriter* writer = new CGNSWriter(solutionFileName.c_str());
            for (int iblock = 0; iblock < blocks.size(); ++iblock)
            {
                const Block* block = blocks[iblock];
                const BlockDef& bdef = blockRankMap.find(block->ID())->second;
                std::cout << "Writing Block " << block->ID() << " solution into CGNS Zone " << bdef.Z << " of " << solutionFileName << std::endl;
                writer->WriteFlowSolution(bdef.Z, *block, block->U(), *Physics::GetInstance());
                writer->WriteTurbulenceSolution(bdef.Z, *block, block->U(), block->UT(), *Physics::GetInstance(), "KOmega");
            }
            delete writer;
        }
        COMM->Barrier();
    }
    double dt = PROF->CheckPoint();
    COMM->Console() << "Saved solution (in cgns format) and it took " << dt << " seconds" << std::endl;
}

void WriteSolution(const char* solutionStem, const char* residualStem, int timeStep,
    const FlowModel& flowModel, const Blocks& blocks, const BlockRankMap& blockRankMap,
    const std::vector<Integrator*>& integrators)
{
    Communicator* COMM = Communicator::GetInstance();
    Profiler* PROF = Profiler::GetInstance();

    std::vector<std::string> fileNames, fileNamesRes;

    int nBlocks = Roster::GetInstance()->GetNumberOfBlocks();

    PROF->CheckPoint();
    for (int rank = 0; rank < COMM->Size(); ++rank)
    {
        for (int iblock = 0; iblock < blocks.size(); ++iblock)
        {
            const Block* block = blocks[iblock];
            std::ostringstream oss;
            oss << solutionStem 
                << '.' << std::setw(std::log10(nBlocks) + 1) << std::setfill('0') << block->ID()
                << '.' << std::setw(5) << timeStep << ".vts";
            fileNames.push_back(oss.str());
            if (rank == COMM->MyRank())
            {
                VTKWriter writer(oss.str().c_str(), VTKWriter::XML);
                writer.AddMesh(block->XYZ());
                writer.AddData(block->U(), 0, "Density", VTKWriter::SCALAR);
                writer.AddData(block->U(), 1, "Momentum", VTKWriter::VECTOR);
                writer.AddData(block->U(), 4, "TotalEnergy", VTKWriter::SCALAR);
                writer.AddData(block->MuK(), 0, "Mu", VTKWriter::SCALAR);
                writer.AddData(block->TurMuK(), 0, "TurMu", VTKWriter::SCALAR);
                writer.AddData(block->UT(), 0, "RhoK", VTKWriter::SCALAR);
                writer.AddData(block->UT(), 1, "RhoOmega", VTKWriter::SCALAR);
                writer.AddData(block->Radius(), 3, "ROmegaSq", VTKWriter::SCALAR);
                writer.Write();
            }
        }

        if (residualStem != NULL)
        {
            for (int iblock = 0; iblock < blocks.size(); ++iblock)
            {
                const Block* block = blocks[iblock];
                std::ostringstream oss;
                oss << residualStem << "." << block->ID() << ".vts";
                if (rank == COMM->MyRank())
                {
                    VTKWriter residWriter(oss.str().c_str(), VTKWriter::XML);
                    Integrator* integ = integrators[iblock];
                    const Structured<double>& resid = integ->ResidualVector();
                    residWriter.AddMesh(block->XYZ());
                    residWriter.AddData(resid, 0, "ResidualDensity", VTKWriter::SCALAR);
                    residWriter.AddData(resid, 1, "ResidualMomentum", VTKWriter::VECTOR);
                    residWriter.AddData(resid, 4, "ResidualTotalEnergy", VTKWriter::SCALAR);
                    residWriter.AddData(integ->DT(), 0, "DT", VTKWriter::SCALAR);
                    residWriter.Write();
                }
            }
        }

        COMM->Barrier();
    }

#if 0
    if (COMM->MyRank() == 0)
    {
        VTKWriter writer("dummy", VTKWriter::XML);
        std::string filename = std::string(solutionStem) + ".vtm";
        writer.Finalize(filename.c_str(), fileNames);
    }
#endif

    double dt = PROF->CheckPoint();
    COMM->Console() << "Saved solution (in vtk format) and it took " << dt << " seconds" << std::endl;
}

void WriteNegativeCells(const char* filename, const std::map<int, std::vector<IndexIJK> >& negativeBlockIndices, double scale)
{
    Communicator* COMM = Communicator::GetInstance();

    // Count the total number of points
    size_t numberOfNegativeCells = 0, totalNumberOfNegativeCells;
    for (std::map<int, std::vector<IndexIJK> >::const_iterator i = negativeBlockIndices.begin();
        i != negativeBlockIndices.end(); ++i)
    {
        numberOfNegativeCells += i->second.size();
    }
    COMM->Console() << "Number of negative cells in Rank " << COMM->MyRank() << " = " << numberOfNegativeCells << std::endl;
    totalNumberOfNegativeCells = COMM->ReduceSum(numberOfNegativeCells, 0);
    if (COMM->MyRank() == 0)
    {
        COMM->Console() << "Total number of negative cells = " << totalNumberOfNegativeCells << std::endl;
    }

    std::ofstream* of;
    COMM->Barrier();
    for (int rank = 0; rank < COMM->Size(); ++rank)
    {
        if (rank == COMM->MyRank())
        {
            if (rank == 0)
            {
                of = new std::ofstream(filename);
                *of << "# vtk DataFile Version 2.0" << std::endl
                    << "Negative Cells" << std::endl
                    << "ASCII" << std::endl
                    << "DATASET POLYDATA" << std::endl
                    << "POINTS " << totalNumberOfNegativeCells << " double" << std::endl;
            }
            else
            {
                of = new std::ofstream(filename, std::ofstream::out | std::ofstream::app);
            }

            for (std::map<int, std::vector<IndexIJK> >::const_iterator i = negativeBlockIndices.begin();
                i != negativeBlockIndices.end(); ++i)
            {
                int blockID = i->first;
                const std::vector<IndexIJK>& indices = i->second;
                // indices contained in "indices" are all from the local blocks, so no need to synchronize here
                Block* block = dynamic_cast<Block*>(Roster::GetInstance()->GetBlock(blockID));
                const Structured<double>& XYZ = block->XYZ();
                for (std::vector<IndexIJK>::const_iterator iv = indices.begin();
                    iv != indices.end(); ++iv)
                {
                    const IndexIJK& ijk = *iv;
                    Vector3 v1(XYZ(ijk - IndexIJK(1, 1, 1)));
                    Vector3 v2(XYZ(ijk - IndexIJK(0, 1, 1)));
                    Vector3 v3(XYZ(ijk - IndexIJK(0, 0, 1)));
                    Vector3 v4(XYZ(ijk - IndexIJK(1, 0, 1)));
                    Vector3 v5(XYZ(ijk - IndexIJK(1, 1, 0)));
                    Vector3 v6(XYZ(ijk - IndexIJK(0, 1, 0)));
                    Vector3 v7(XYZ(ijk - IndexIJK(0, 0, 0)));
                    Vector3 v8(XYZ(ijk - IndexIJK(1, 0, 0)));
                    Vector3 vc = 0.125 * (v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8) / scale;
                    *of << vc.X() << ' ' << vc.Y() << ' ' << vc.Z() << std::endl;
                }
            }
            of->close();
            delete of;
        }
        COMM->Barrier();
    }
}

int main(int argc, char** argv)
{
    Communicator::Initialize(&argc, &argv);
    Communicator* COMM = Communicator::GetInstance();
    int myrank = COMM->MyRank();
    std::ostream& CONSOLE = COMM->Console();

    if (myrank == 0)
    {
        splash(std::cout);
    }

    assert(argc > 1);

    bool probe = false;
    BlockRanges probeRanges;
    if (argc > 2)
    {
        CONSOLE << "Reading block ranges for debugging output from file " << argv[2] << std::endl;
        ReadBlockRanges(probeRanges, argv[2]);
        probe = probeRanges.size() > 0;
    }

    Blocks blocks; // local blocks
    double meshScale = 1.0;
    bool unsteady = false;
    double dTimeReal, Time0;
    int subIterations = 1;
    int orderMUSCL;
    int solverType = 0; // LUSGS
    double cfl0 = 2.0, cflmax = 1000.0, ser = 1.0;
    bool localTimeStepping;
    int maxIter, iterT0, iterStartSER = 5;
    BlockRankMap blockRankMap;
    typedef std::vector<InterfaceDataExchanger*> InterfaceDataExchangers;
    InterfaceDataExchangers IDXsU, IDXsUT;
    ResidualEvaluatorBase* resEval = NULL;
    FlowModel* model;
    KOmega1988::TurbulenceModel* modelT;
    KOmega1988::ResidualEvaluator* resEvalT = NULL;

    CGNSStructure cgnsSt;

    for (int rank = 0; rank < COMM->Size(); ++rank)
    {
        if (rank != myrank)
        {
            COMM->Barrier();
            continue;
        }

        int dummy;

        std::ifstream f(argv[1]);
        std::string line;
        std::istringstream iss;

        // Mesh file
        std::getline(f, line); std::getline(f, line);
        std::string meshFileName = line;
        CGNSReader meshReader(meshFileName.c_str());
        cgnsSt = meshReader.ReadStructure();

        // Scale
        std::getline(f, line); std::getline(f, line);
        iss.clear(); iss.str(line);
        iss >> meshScale;

        // Reference values
        std::getline(f, line); std::getline(f, line);
        iss.clear(); iss.str(line);
        double gamma, rhoRef, TRef, RGAS;
        iss >> gamma >> rhoRef >> TRef >> RGAS;
        Physics::Initialize(gamma, rhoRef, TRef, RGAS);
        std::cout << "Reynolds Number (based on reference acoustic speed) = " << Physics::GetInstance()->ReynoldsNumber() << std::endl;

        // Time advancement
        std::getline(f, line); std::getline(f, line);
        iss.clear(); iss.str(line);
        iss >> dummy;
        unsteady = dummy != 0;
    
        std::getline(f, line);
        iss.clear(); iss.str(line);
        if (unsteady)
        {
            iss >> dTimeReal >> Time0;
            std::cout << "Unsteady simulation, dTimeReal = " << dTimeReal << ", Time0 = " << Time0 << std::endl;
        }

        // Solver
        std::getline(f, line); std::getline(f, line);
        iss.clear(); iss.str(line);
        int local = 1;
        iss >> solverType >> cfl0 >> cflmax >> ser >> iterStartSER >> local;
        localTimeStepping = local != 0;

        // Convective flux
        std::getline(f, line); std::getline(f, line);
        iss.clear(); iss.str(line);
        std::string convectiveFluxType, limiterType;
        iss >> convectiveFluxType >> orderMUSCL >> limiterType;
        if (orderMUSCL == 1)
        {
            std::cout << "First order Roe flux" << std::endl;
            CONSOLE << "First order Roe flux" << std::endl;
            resEval = new ::ResidualEvaluator<5, FirstOrderReconstructor<5> >();
            resEvalT = new KOmega1988::ResidualEvaluator();
        }
        else if (orderMUSCL > 1)
        {
            std::cout << "Second order Roe flux with minmod limiter" << std::endl;
            CONSOLE << "Second order Roe flux with minmod limiter" << std::endl;
            resEval = new ::ResidualEvaluator<5, Reconstructor<5, MinModLimiter, PrimitiveVariableCodec> >();
            resEvalT = new KOmega1988::ResidualEvaluator(); // FIXME: 1st-order
        }
        else
        {
            CONSOLE << "orderMUSCL = " << orderMUSCL << " not supported." << std::endl;
            throw 666;
        }
        model = new FlowModel;
        modelT = new KOmega1988::TurbulenceModel();

        // Iterations
        std::getline(f, line); std::getline(f, line);
        iss.clear(); iss.str(line);
        iss >> maxIter >> iterT0 >> subIterations;
        if (!unsteady)
        {
            subIterations = 1;
        }
        std::cout << maxIter << " iterations" << std::endl;
        std::cout << "At iteration " << iterT0 << " turbulence equations will be activated." << std::endl;
        std::cout << "Unsteady/steady = " << unsteady << ", subIterations = " << subIterations << std::endl;

        // # of blocks
        std::getline(f, line); std::getline(f, line);
        iss.clear(); iss.str(line);
        int numBlocks;
        iss >> numBlocks;

        // Blocks. We create blocks here already.
        std::getline(f, line);
        for (int i = 0; i < numBlocks; ++i)
        {
            // Block dimension
            int blockID, rank, Z, imin, jmin, kmin, imax, jmax, kmax;
            std::getline(f, line); // comment line
            std::getline(f, line);
            iss.clear(); iss.str(line);
            iss >> blockID >> rank >> Z >> imin >> jmin >> kmin >> imax >> jmax >> kmax;
            // Rigid body motion
            std::string rbmType;
            RigidBodyMotion* rbm = NULL;
            BlockOffset offset;
            std::getline(f, line); // comment line
            std::getline(f, line);
            iss.clear(); iss.str(line);
            iss >> rbmType;
            if (rbmType == "ROTATIONAL")
            {
                double ox, oy, oz;
                iss >> ox >> oy >> oz;
                rbm = new RotationalMotion(Vector3(0.0, 0.0, 0.0), Vector3(ox, oy, oz));
                double oix, oiy, oiz; // initial offset
                iss >> oix >> oiy >> oiz;
                offset = BlockOffset(BlockOffset::ROTATION, Vector3(oix, oiy, oiz));
            }
            else if (rbmType == "TRANSLATIONAL")
            {
                assert(false);
            }
            // Periodicity
            std::string perType;
            Vector3 periodicity;
            std::getline(f, line); // comment line
            std::getline(f, line);
            iss.clear(); iss.str(line);
            iss >> perType;
            if (perType == "ROTATIONAL")
            {
                double px, py, pz;
                iss >> px >> py >> pz;
                periodicity = Vector3(px, py, pz);
            }
            else if (perType == "TRANSLATIONAL")
            {
                assert(false);
            }

            BlockDef bdef(blockID, rank, Z, IndexRange(imin, jmin, kmin, imax, jmax, kmax), offset);
            VirtualBlock* vblock;
            if (bdef.Rank == myrank)
            {
                vblock = Block::New(bdef.BlockID, bdef.MeshRange, unsteady);
                Block* block = dynamic_cast<Block*>(vblock);
                bdef.LocalBlock = block;
                blocks.push_back(block);
            }
            else
            {
                vblock = VirtualBlock::New(bdef.BlockID, bdef.Rank, bdef.MeshRange);
            }
            blockRankMap[blockID] = bdef;
            if (rbm != NULL)
            {
                vblock->SetRigidBodyMotion(*rbm);

                // Record this info in a CGNSStructure so that solution output files can contain this info later.
                CGNSStructure::Zone& cgnsZone = cgnsSt.Bases[0].FindZone(Z);
                cgnsZone.Motion = rbm->Clone();

                delete rbm;
            }
            vblock->SetPeriodicity(periodicity);
        }

        // Block patches
        std::getline(f, line); std::getline(f, line);
        iss.clear(); iss.str(line);
        int numBlockPatches;
        iss >> numBlockPatches;
        for (int i = 0; i < numBlockPatches; ++i)
        {
            int blockID, is, js, ks, ie, je, ke, uniqueID;
            std::string blockPatchName;
            std::getline(f, line);
            iss.clear(); iss.str(line);
            iss >> blockID >> is >> js >> ks >> ie >> je >> ke >> uniqueID >> blockPatchName;
            IndexRange mr(is, js, ks, ie, je, ke);
            BlockPatch::New(blockID, mr, uniqueID);
        }

        // Block patch families
        std::getline(f, line); std::getline(f, line);
        iss.clear(); iss.str(line);
        int numBlockPatchFamilies;
        iss >> numBlockPatchFamilies;
        for (int i = 0; i < numBlockPatchFamilies; ++i)
        {
            int familyID, blockPatchID;;
            BlockPatches blockPatches;
            std::getline(f, line);
            iss.clear(); iss.str(line);
            iss >> familyID;
            while (iss >> blockPatchID)
                blockPatches.push_back(Roster::GetInstance()->GetBlockPatch(blockPatchID));
            Roster::GetInstance()->RegisterBlockPatchFamily(familyID, blockPatches);
        }
        const Roster::BlockPatchesMap& bpFams = Roster::GetInstance()->BlockPatchFamilies();
        for (Roster::BlockPatchesMap::const_iterator i = bpFams.begin(); i != bpFams.end(); ++i)
        {
            int fam = i->first;
            const BlockPatches& bps = i->second;
            CONSOLE << "Patch Family " << fam << std::endl;
            for (BlockPatches::const_iterator j = bps.begin(); j != bps.end(); ++j)
            {
                CONSOLE << "    " << *j << std::endl;
            }
        }

        // Read mesh
        for (BlockRankMap::iterator i = blockRankMap.begin(); i != blockRankMap.end(); ++i)
        {
            BlockDef& bdef = i->second;

            if (bdef.Rank == myrank)
            {
                std::cout << "Reading Block " << bdef.BlockID << ": Zone = " << bdef.Z
                    << ", MeshRange = " << bdef.MeshRange << " " << bdef.MeshRange.Count() << " vertices" << std::endl;
                CONSOLE << "Reading Block " << bdef.BlockID << ": Zone = " << bdef.Z
                    << ", MeshRange = " << bdef.MeshRange << " " << bdef.MeshRange.Count() << " vertices" << std::endl;

                Block* block = dynamic_cast<Block*>(Roster::GetInstance()->GetBlock(bdef.BlockID));
                assert(block != NULL);
                std::string zoneName;
                meshReader.ReadMesh(bdef.Z, block->XYZ(), zoneName, bdef.MeshRange);
                CONSOLE << "read zone " << bdef.Z << " named " << zoneName << std::endl;
                CONSOLE << "scaling the mesh by " << meshScale << std::endl;
                block->XYZ() *= meshScale;
                CONSOLE << "offseting the mesh" << std::endl;
                bdef.Offset.ApplyOffset(block->XYZ()); // NOTE: offset is applied *after* scaling
                block->ComputeMetrics();
            }
        }

        // Local zone
        for (BlockRankMap::const_iterator i = blockRankMap.begin(); i != blockRankMap.end(); ++i)
        {
            int Z = i->first;
            const BlockDef& bdef = i->second;
            bool local = bdef.Rank == myrank;

            std::cout << "**************" << std::endl;
            std::cout << "Setting up Block " << bdef.BlockID << " for Rank " << bdef.Rank << std::endl;
            std::cout << "**************" << std::endl;
            if (!local)
            {
                CONSOLE << "Skipping Block " << bdef.BlockID << " since it's not local." << std::endl;
                std::getline(f, line); // # Zone
                std::getline(f, line);

                std::getline(f, line); // # BCs
                std::getline(f, line);
                int numBCs;
                iss.clear(); iss.str(line); iss >> numBCs;
                for (int ibc = 0; ibc < numBCs; ++ibc)
                {
                    std::getline(f, line);
                }

                std::getline(f, line); // # 1-to-1 connectivities
                std::getline(f, line);
                int numConn1to1s;
                iss.clear(); iss.str(line); iss >> numConn1to1s;
                for (int iconn = 0; iconn < numConn1to1s; ++iconn)
                {
                    std::getline(f, line);
                }

                std::getline(f, line); // # Non-matching connectivities
                std::getline(f, line);
                int numConnNMBs;
                iss.clear(); iss.str(line); iss >> numConnNMBs;
                for (int iconn = 0; iconn < numConnNMBs; ++iconn)
                {
                    std::getline(f, line);
                }

                std::getline(f, line); // # Interfaces
                std::getline(f, line);
                int numIFs;
                iss.clear(); iss.str(line); iss >> numIFs;
                for (int iif = 0; iif < numIFs; ++iif)
                {
                    std::getline(f, line);
                    for (int isd = 0; isd < 2; ++isd)
                    {
                        std::getline(f, line);
                        int npatches;
                        iss.clear(); iss.str(line); iss >> npatches;
                        for (int np = 0; np < npatches; ++np)
                        {
                            std::getline(f, line);
                        }
                    }
                }

                std::getline(f, line); // # Flowfield Initialization
                std::getline(f, line);
                int numSpecs;
                std::getline(f, line);
                iss.clear(); iss.str(line); iss >> numSpecs;
                for (int iinit = 0; iinit < numSpecs; ++iinit)
                {
                    std::getline(f, line);
                }
                std::getline(f, line); // # Turbulence Initialization
                std::getline(f, line);
                std::getline(f, line);
                continue;
            }

            CONSOLE << "**************" << std::endl
                << "Setting up Block " << bdef.BlockID << " for Rank " << bdef.Rank << std::endl
                << "**************" << std::endl;

            std::getline(f, line); // # Zone
            CONSOLE << '"' << line << '"' << std::endl;
            std::getline(f, line);
            CONSOLE << '"' << line << '"' << std::endl;
            iss.clear(); iss.str(line);
            iss >> Z;
            CONSOLE << "Z (was just read) = " << Z << " bdef.Z = " << bdef.Z << std::endl;
            Block* block = bdef.LocalBlock;

            // BCs
            int numBCs;
            std::getline(f, line);
            std::getline(f, line);
            iss.clear(); iss.str(line); iss >> numBCs;

            std::cout << "Reading BCs for Zone " << Z << std::endl;
            CONSOLE << "Reading BCs for Zone " << Z << std::endl;
            for (int ibc = 0; ibc < numBCs; ++ibc)
            {
                std::string name, type;
                std::string strData;
                double data[5];
                TurbulenceSpec* turbSpec = NULL;
                std::getline(f, line);

                size_t i1, i2;
                i1 = line.find_first_of('\"', 0);
                i2 = line.find_first_of('\"', i1 + 1);
                name = line.substr(i1 + 1, i2 - i1 - 1);
                line = &line[i2 + 1];

                iss.clear(); iss.str(line);

                // Family ID and Patch ID
                int familyID, bpID;
                iss >> familyID >> bpID;

                // Range
                int r[6];
                iss >> r[0] >> r[1] >> r[2] >> r[3] >> r[4] >> r[5];
                IndexRange range(r[0], r[1], r[2], r[3], r[4], r[5]);
                const BlockPatch& bcbp = Roster::GetInstance()->GetBlockPatch(bdef.BlockID, range);
                CONSOLE << "DEBUG: " << name << " " << range << " is BlockPatch " << bcbp.UniqueID() << " " << bcbp.MeshRange() << std::endl;

                // Type
                iss >> type;
                if (type == "Fix")
                {
                    iss >> data[0] >> data[1] >> data[2] >> data[3] >> data[4];
                    turbSpec = ReadTurbulenceSpec(iss);
                }
                else if (type == "OutletStaticPressure")
                {
                    iss >> data[0];
                }
                else if (type == "ViscousWall")
                {
                    iss >> strData >> data[0];
                }
                else if (type == "InletTotal")
                {
                    iss >> data[0] >> data[1] >> data[2] >> data[3] >> data[4];
                    turbSpec = ReadTurbulenceSpec(iss);
                }
                std::cout << "  BC Name = " << name << ", range = " << range << std::endl;
                CONSOLE << "  BC Name = " << name << ", range = " << range << std::endl;

                Direction dir = IndexUtils::PatchDirection(block->MeshRange(), range);
                BC* bc;
                if (type == "Symmetry" || type == "InviscidWall")
                {
                    bc = new BCSymmetry(range, dir);
                }
                else if (type == "Fix")
                {
                    double* U0 = &data[0];
                    bc = new BCFix(range, dir, block->U().DOF(), U0, 2, turbSpec);
                }
                else if (type == "Extrapolate")
                {
                    bc = new BCExtrapolate(range, dir);
                }
                else if (type == "ViscousWall")
                {
                    bool isStationary = strData == "Stationary";
                    bc = new BCViscousWall(range, dir, isStationary, data[0]);
                }
                else if (type == "RotatingViscousWall")
                {
                    bc = new BCRotatingViscousWall(range, dir, data[0], Vector3(data[1], data[2], data[3]));
                }
                else if (type == "OutletStaticPressure")
                {
                    bc = new BCOutletStaticPressure(range, dir, data[0]);
                }
                else if (type == "InletTotal")
                {
                    bc = new BCInletTotal(range, dir, data[0], data[1], Vector3(data[2], data[3], data[4]), 2, turbSpec);
                }
                else if (type == "ForcePull")
                {
                    bc = new BCForcePull(*model, range, dir);
                }
                else
                {
                    assert(false);
                }
                block->RegisterBC(bc);
                std::cout << "Added BC " << name << " of type " << type << bc->MeshRange() << " : " << range << std::endl;
                CONSOLE << "Added BC " << name << " of type " << type << bc->MeshRange() << " : " << range << std::endl;

                delete turbSpec; // the TurbulenceSpec object passed to the BC ctors above have been Clone'd by them, so no need to keep this around.
            }

            // Connectivities
            int numConn1to1s;
            std::getline(f, line);
            std::getline(f, line);
            iss.clear(); iss.str(line); iss >> numConn1to1s;

            std::cout << "Reading 1-to-1 connectivities for Block " << Z << std::endl;
            CONSOLE << "Reading 1-to-1 connectivities for Block " << Z << std::endl;
            for (int iconn = 0; iconn < numConn1to1s; ++iconn)
            {
                std::getline(f, line);
                iss.clear(); iss.str(line);

                std::string connName;
                int tag;
                int donorBlock;
                iss >> connName >> tag >> donorBlock;

                int sr[6], dr[6], tr[3];
                double center[3], angle[3], translation[3];
                iss >> sr[0] >> sr[1] >> sr[2] >> sr[3] >> sr[4] >> sr[5];
                iss >> dr[0] >> dr[1] >> dr[2] >> dr[3] >> dr[4] >> dr[5];
                iss >> tr[0] >> tr[1] >> tr[2]
                    >> center[0] >> center[1] >> center[2]
                    >> angle[0] >> angle[1] >> angle[2]
                    >> translation[0] >> translation[1] >> translation[2];

                IndexRange selfRange(sr[0], sr[1], sr[2], sr[3], sr[4], sr[5]);
                IndexRange donorRange(dr[0], dr[1], dr[2], dr[3], dr[4], dr[5]);
                IndexTransform::Canonize(tr, selfRange, donorRange);
                IndexTransform transform(tr, selfRange, donorRange);

                Direction dir = IndexUtils::PatchDirection(block->MeshRange(), selfRange);
                Connectivity1to1* conn = new Connectivity1to1(
                    selfRange, dir, *block, transform, donorBlock, tag
                    );

                /// Periodicity
                bool angularPerio = (angle[0] != 0.0 || angle[1] != 0.0 || angle[2] != 0.0);
                bool translationalPerio = (translation[0] != 0.0 || translation[1] != 0.0 || translation[2] != 0.0);
                if (angularPerio)
                {
                    Vector3 angleInRad(angle);
                    //angleInRad *= M_PI / 180.0;
                    conn->SetPeriodicity(Connectivity1to1::ROTATION, Vector3(center), angleInRad);
                }
                else if (translationalPerio)
                {
                    //assert(false);
                }
                block->RegisterBC(conn);
                std::cout << "Added Conn1to1 " << connName << " : self " << selfRange << " donor " << donorRange << std::endl;
                CONSOLE << "Added Conn1to1 " << connName << " : self " << selfRange << " donor " << donorRange << std::endl;
            }

            int numConnNMBs;
            std::getline(f, line);
            std::getline(f, line);
            iss.clear(); iss.str(line); iss >> numConnNMBs;

            std::cout << "Reading non-matching connectivities for Block " << Z << std::endl;
            CONSOLE << "Reading non-matching connectivities for Block " << Z << std::endl;
            for (int iconn = 0; iconn < numConnNMBs; ++iconn)
            {
                std::getline(f, line);
                iss.clear(); iss.str(line);

                std::string connName;
                int tag;
                int donorBlock;
                iss >> connName >> tag >> donorBlock;

                int sr[6], dr[6], tr[3];
                double center[3], angle[3], translation[3];
                iss >> sr[0] >> sr[1] >> sr[2] >> sr[3] >> sr[4] >> sr[5];
                iss >> dr[0] >> dr[1] >> dr[2] >> dr[3] >> dr[4] >> dr[5];
                iss >> tr[0] >> tr[1] >> tr[2];

                IndexRange selfRange(sr[0], sr[1], sr[2], sr[3], sr[4], sr[5]);
                IndexRange donorRange(dr[0], dr[1], dr[2], dr[3], dr[4], dr[5]);
                IndexTransform::Canonize(tr, selfRange, donorRange);
                IndexTransform transform(tr, selfRange, donorRange);

                Direction dir = IndexUtils::PatchDirection(block->MeshRange(), selfRange);
                ConnectivityAbutting* conn = new ConnectivityAbutting(
                    selfRange, dir, *block, donorRange, donorBlock, tag
                    );
            }

            // Interfaces
            int numIFs;
            std::getline(f, line);
            std::getline(f, line);
            iss.clear(); iss.str(line); iss >> numIFs;
            std::cout << "Reading interfaces for Block " << Z << std::endl;
            CONSOLE << "Reading interfaces for Block " << Z << std::endl;
            for (int iif = 0; iif < numIFs; ++iif)
            {
                std::string ifname;
                std::getline(f, line);
                iss.clear(); iss.str(line); iss >> ifname;

                int numSelfPatches, numDonorPatches;
                BlockPatches selfPatches, donorPatches;
                std::getline(f, line);
                iss.clear(); iss.str(line); iss >> numSelfPatches;
                std::cout << "there are " << numSelfPatches << " self patches" << std::endl;
                CONSOLE << "there are " << numSelfPatches << " self patches" << std::endl;
                for (int ip = 0; ip < numSelfPatches; ++ip)
                {
                    std::string dummy, patchName;
                    int blockID, r[6];
                    std::getline(f, line);
                    iss.clear(); iss.str(line);
                    iss >> dummy >> blockID >> r[0] >> r[1] >> r[2] >> r[3] >> r[4] >> r[5];
                    patchName = UnQuote(dummy);
                    IndexRange range(r[0], r[1], r[2], r[3], r[4], r[5]);
                    const BlockPatch& bp = Roster::GetInstance()->GetBlockPatch(blockID, range);
                    selfPatches.push_back(bp);
                }
                std::getline(f, line);
                iss.clear(); iss.str(line); iss >> numDonorPatches;
                std::cout << "there are " << numDonorPatches << " donor patches" << std::endl;
                CONSOLE << "there are " << numDonorPatches << " donor patches" << std::endl;
                for (int ip = 0; ip < numDonorPatches; ++ip)
                {
                    std::string dummy, patchName;
                    int blockID, r[6];
                    std::getline(f, line);
                    iss.clear(); iss.str(line);
                    iss >> dummy >> blockID >> r[0] >> r[1] >> r[2] >> r[3] >> r[4] >> r[5];
                    patchName = UnQuote(dummy);
                    IndexRange range(r[0], r[1], r[2], r[3], r[4], r[5]);
                    const BlockPatch& bp = Roster::GetInstance()->GetBlockPatch(blockID, range);
                    donorPatches.push_back(bp);
                }
                std::cout << "read self-donor = " << selfPatches.size() << ", " << donorPatches.size() << std::endl;
                CONSOLE << "read self-donor = " << selfPatches.size() << ", " << donorPatches.size() << std::endl;
                AbuttingInterface* interface = SimplePlanarAbuttingInterface::New(selfPatches, donorPatches);

                // Read the mesh for blocks affected by the interface.
                std::set<int> blockIDs = interface->BlockIDs();
                for (std::set<int>::const_iterator i = blockIDs.begin(); i != blockIDs.end(); ++i)
                {
                    int blockID = *i;
                    CONSOLE << "Reading block " << blockID << " mesh (for abutting interface " << ifname << ")" << std::endl;
                    const BlockDef& bdefDonor = blockRankMap[blockID];
                    Structured<double> XYZ(3, bdefDonor.MeshRange);
                    std::string zoneName;
                    meshReader.ReadMesh(bdefDonor.Z, XYZ, zoneName);
                    CONSOLE << "scaling the mesh by " << meshScale << std::endl;
                    XYZ *= meshScale;
                    CONSOLE << "offseting the mesh" << std::endl;
#if 0
                    std::ostringstream fileNameSS;
                    fileNameSS << "mesh." << bdef.BlockID << "." << bdefDonor.BlockID << ".before.vtk";
                    VTKWriter* meshWriter = new VTKWriter(fileNameSS.str().c_str());
                    meshWriter->AddMesh(XYZ);
                    meshWriter->Write();
                    delete meshWriter;
#endif
                    bdefDonor.Offset.ApplyOffset(XYZ); // NOTE: offset is applied *after* scaling
#if 0
                    fileNameSS.str("");
                    fileNameSS << "mesh." << bdef.BlockID << "." << bdefDonor.BlockID << ".after.vtk";
                    meshWriter = new VTKWriter(fileNameSS.str().c_str());
                    meshWriter->AddMesh(XYZ);
                    meshWriter->Write();
                    delete meshWriter;
#endif
                    interface->SetPatchMesh(blockID, XYZ);
                    delete[] XYZ.Data;
                }

                // Create data exchangers (each for U and UT) for this interface
                IDXsU.push_back(new InterfaceDataExchanger(*interface, model, "U"));
                IDXsUT.push_back(new InterfaceDataExchanger(*interface, modelT, "UT"));
            }

            // Initialization
            Physics* PHYS = Physics::GetInstance();

            std::string initMethod = "UNIFORM";
            int numSpecs;

            std::getline(f, line);
            std::getline(f, line);
            iss.clear(); iss.str(line);
            iss >> initMethod;
            std::getline(f, line);
            iss.clear(); iss.str(line);
            iss >> numSpecs;

            if (initMethod == "UNIFORM")
            {
                assert(numSpecs == 1);
                double U0[5], U0Local[5];
                std::getline(f, line);
                iss.clear(); iss.str(line);
                iss >> U0[0] >> U0[1] >> U0[2] >> U0[3] >> U0[4];
                U0[0] /= PHYS->RhoRef();
                U0[1] /= PHYS->RhoRef() * PHYS->VRef();
                U0[2] /= PHYS->RhoRef() * PHYS->VRef();
                U0[3] /= PHYS->RhoRef() * PHYS->VRef();
                U0[4] /= PHYS->RhoRef() * PHYS->ERef();
                Initialize(*model, *block, block->U(), U0); //block->U().SetTo(U0);
            }
            else if (initMethod == "ISENTROPICEXPANSION")
            {
                std::string dir;
                double P0, T0, pos, pres;
                std::getline(f, line);
                iss.clear(); iss.str(line);
                iss >> dir >> P0 >> T0;
                IsentropicExpansionInitializer initer(*model, dir.c_str());
                initer.SetTotalQuantities(P0, T0);
                for (int i = 0; i < numSpecs - 1; ++i)
                {
                    std::getline(f, line);
                    iss.clear(); iss.str(line);
                    iss >> pos >> pres;
                    CONSOLE << "IsentropicExpansionIniter: pos = " << pos << ", p = " << pres << std::endl;
                    initer.AddPressureSpec(pos, pres);
                }
                initer.Initialize(block->U(), *block);
            }
            else if (initMethod == "CENTRIFUGAL")
            {
                std::string axis;
                double P0, T0, axial, radial, vx, vy, vz, pres;
                std::getline(f, line);
                iss.clear(); iss.str(line);
                iss >> axis >> P0 >> T0;
                CentrifugalInitializer initer(*model, axis.c_str());
                initer.SetTotalQuantities(P0, T0);
                for (int i = 0; i < numSpecs - 1; ++i)
                {
                    std::getline(f, line);
                    iss.clear(); iss.str(line);
                    iss >> axial >> radial >> vx >> vy >> vz >> pres;
                    Vector3 vdir(vx, vy, vz);
                    CONSOLE << "CentrifugalIniter: a = " << axial << ", r = " << radial << ", vdir = " << vdir << ", p = " << pres << std::endl;
                    initer.AddSpec(axial, radial, vdir, pres);
                }
                initer.Initialize(block->U(), *block);
            }
            else if (initMethod == "READFROMFILE")
            {
                assert(numSpecs == 1);
                std::string solFileName;
                std::getline(f, line);
                iss.clear(); iss.str(line);
                iss >> solFileName;
                std::cout << "Reading the initial flow solution from " << solFileName << std::endl;
                CONSOLE << "Reading the initial flow solution from " << solFileName << std::endl;
                CGNSReader solReader(solFileName.c_str());
                solReader.ReadFlowSolution(bdef.Z, *block, block->U(), *Physics::GetInstance());
            }
            else
            {
                CONSOLE << "InitMethod " << initMethod << " not supported" << std::endl;
                throw 666;
            }

            // Turbulence field initialization
            TurbulenceSpec* turbSpec = NULL;
            std::string initMethodT = "UNIFORM";
            std::getline(f, line);
            std::getline(f, line);
            iss.clear(); iss.str(line);
            iss >> initMethodT;

            if (initMethodT == "UNIFORM")
            {
                // FIXME: rewrite this whole block
                double UT0[2];
                std::getline(f, line);
                iss.clear(); iss.str(line);
                turbSpec = ReadTurbulenceSpec(iss);
                CONSOLE << "initMethodT == UNIFORM, TurbSpec = " << *turbSpec << std::endl;
                modelT->SetUT(block->UT(), block->U(), block->CellRange(), *turbSpec);
            }
            else if (initMethodT == "READFROMFILE")
            {
                std::string solFileName;
                std::getline(f, line);
                iss.clear(); iss.str(line);
                iss >> solFileName;
                std::cout << "Reading the initial turbulence solution from " << solFileName << std::endl;
                CONSOLE << "Reading the initial turbulence solution from " << solFileName << std::endl;
                CGNSReader solReader(solFileName.c_str());
                solReader.ReadTurbulenceSolution(bdef.Z, *block, block->U(), block->UT(), *Physics::GetInstance(), "KOmega"); // FIXME: make sure this reads correct turbulence quantities and UT is correctly initialized accordingly.
            }
            else
            {
                CONSOLE << "InitMethodT " << initMethodT << " not supported" << std::endl;
                throw 666;
            }
            delete turbSpec;

            block->CheckNegatives();
            block->ApplyBCs();
            block->ApplyTurbBCs();
        }
        f.close();
        COMM->Barrier();
    }

    // Exchange rind data among blocks
    {
    CONSOLE << "Exchanging and filling rind cell data" << std::endl;
    COMM->Barrier();
    typedef std::vector<StructuredDataExchanger*> SDXs;
    SDXs sdxUs, sdxUTs;
    for (Blocks::iterator i = blocks.begin(); i != blocks.end(); ++i)
    {
        Block* block = *i;
        sdxUs.push_back(new StructuredDataExchanger(*block, block->U()));
        sdxUTs.push_back(new StructuredDataExchanger(*block, block->UT()));
    }
    for (SDXs::iterator i = sdxUs.begin(); i != sdxUs.end(); ++i)
    {
        (*i)->Start();
    }
    COMM->Barrier();
    for (SDXs::iterator i = sdxUs.begin(); i != sdxUs.end(); ++i)
    {
        (*i)->Finish();
    }
    COMM->Barrier();
    for (SDXs::iterator i = sdxUTs.begin(); i != sdxUTs.end(); ++i)
    {
        (*i)->Start();
    }
    COMM->Barrier();
    for (SDXs::iterator i = sdxUTs.begin(); i != sdxUTs.end(); ++i)
    {
        (*i)->Finish();
    }
    COMM->Barrier();
    for (SDXs::iterator i = sdxUs.begin(); i != sdxUs.end(); ++i)
    {
        delete (*i);
    }
    for (SDXs::iterator i = sdxUTs.begin(); i != sdxUTs.end(); ++i)
    {
        delete (*i);
    }
    }

    // Create integrators
    CONSOLE << "Creating integrators" << std::endl;
    Integrators integrators, integratorsT;

    for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
    {
        Block* block = blocks[iblock];
        if (solverType == 0)
        {
            integrators.push_back(new LUSGSIntegrator<FlowModel, ResidualEvaluatorBase>(*model, resEval, *block, block->U(), block->UStorage(1), block->UStorage(2)));
            integratorsT.push_back(new LUSGSIntegrator<KOmega1988::TurbulenceModel, KOmega1988::ResidualEvaluator>(*modelT, resEvalT, *block, block->UT(), block->UTStorage(1), block->UTStorage(2)));
        }
        else if (solverType == 1)
        {
            integrators.push_back(new SimpleExplicitIntegrator<FlowModel, ResidualEvaluatorBase>(*model, resEval, *block, block->U(), block->UStorage(1), block->UStorage(2)));
            integratorsT.push_back(new SimpleExplicitIntegrator<KOmega1988::TurbulenceModel, KOmega1988::ResidualEvaluator>(*modelT, resEvalT, *block, block->UT(), block->UTStorage(1), block->UTStorage(2)));
        }
    }

    double cfl = cfl0;
    TimeStepEvaluator<FlowModel> tsEval(*model, cfl, localTimeStepping);
    TimeStepEvaluator<KOmega1988::TurbulenceModel> tsEvalT(*modelT, cfl0, localTimeStepping);

    COMM->Barrier();

    std::ofstream frms("rms.dat");

    Residual res0, res0T;
    double rms0[5], rmsTotal0, rms0T[2], rmsTotal0T;

    //Probe(blocks);

    IterationContext iteration;
    if (unsteady)
    {
        iteration.SetType(IterationContext::UNSTEADY);
        iteration.SetDT(dTimeReal);
        iteration.AdvanceTimeStep(); // the very first solution we are solving for is at t = t0 + dt
        cfl0 = cflmax;
    }

    // Initialize interface mapping
    CONSOLE << "Initializing interface mappers" << std::endl;
    Interfaces& interfaces = Roster::GetInstance()->GetInterfaces();
    for (Interfaces::iterator i = interfaces.begin(); i != interfaces.end(); ++i)
    {
        (*i)->InitializeMapping();
        (*i)->MapMesh(iteration);
    }

    // Exchange data across interfaces
    CONSOLE << "Exchanging data across interfaces" << std::endl;
    for (InterfaceDataExchangers::iterator iidx = IDXsU.begin(); iidx != IDXsU.end(); ++iidx)
    {
        (*iidx)->Exchange();
    }
    COMM->Barrier();

    for (Blocks::iterator i = blocks.begin(); i != blocks.end(); ++i)
    {
        Block* block = *i;

        // 2D flow enforcement
        //block->RegisterBC(new BCEnforce2D(BCEnforce2D::TWOD_Z));

        block->ApplyBCs();
        block->ComputeTransportProperties();
        block->ApplyTurbBCs();
        block->ComputeTurbulentTransportProperties();
        block->CheckNegatives();
    }

    if (unsteady)
    {
        // FIXME: (to verify) We need to initialize two(?) extra solution layers.
        for (Blocks::iterator i = blocks.begin(); i != blocks.end(); ++i)
        {
            Block* block = *i;
            //block->ApplyBCs();
            //block->ComputeTransportProperties();
            block->ShiftTime();
            block->ShiftTime();
        }
    }

    Probe(CONSOLE, blocks, probeRanges);

    CONSOLE << "STARTING THE OUTER LOOP" << std::endl;
    for (int iloop = 0; iloop < maxIter; ++iloop)
    {
        if (COMM->MyRank() == 0)
        {
            std::cout << "Iteration " << iloop << '\t' << cfl << std::endl;
        }

        bool bailout = false;

        if (unsteady)
        {
            // Map the mesh interface
            Interfaces& interfaces = Roster::GetInstance()->GetInterfaces();
            for (Interfaces::iterator i = interfaces.begin(); i != interfaces.end(); ++i)
            {
                (*i)->MapMesh(iteration);
            }
        }

        for (int isub = 0; isub < subIterations; ++isub)
        {
            // ******************************************
            // Flow Solve
            // ******************************************
            tsEval.Evaluate(integrators);

            int nisteps = integrators[0]->IntegrationSteps();
            for (int ii = 0; ii < nisteps; ++ii)
            {
                // Pre-Integrate Phase
                for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                {
                    integrators[iblock]->PreIntegrateStart(ii, iteration);
                }
                for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                {
                    integrators[iblock]->PreIntegrateFinish(ii, iteration);
                }

                // Integrate Phase
                for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                {
                    Integrator* integrator = integrators[iblock];
                    integrator->Integrate(ii, iteration, blocks[iblock]->U());
                }
                if (ii == nisteps - 1)
                {
                    for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                    {
                        Block* block = blocks[iblock];
                        Integrator* integrator = integrators[iblock];
                        //block->U().Add(1.0, integrator->DU(), block->CellRange());
                        model->UplusDU(block->U(), integrator->DU(), block->CellRange(), *block);
                        block->ComputeTransportProperties();
                    }
                }

                // Post-Integrate Phase
                for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                {
                    integrators[iblock]->PostIntegrateStart(ii, iteration);
                }
                for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                {
                    integrators[iblock]->PostIntegrateFinish(ii, iteration);
                }
                COMM->Barrier();
            }

            // Interface data exchange
            for (InterfaceDataExchangers::iterator iidx = IDXsU.begin(); iidx != IDXsU.end(); ++iidx)
            {
                (*iidx)->Exchange();
            }

            // Accumulate the residuals from the local blocks
            Residual res = integrators[0]->LatestResidual();
            for (size_t iblock = 1; iblock < blocks.size(); ++iblock)
            {
                res += integrators[iblock]->LatestResidual();
            }

            // Local residual statistics
            double rms[5], rmsTotal;
            res.GetRMS(rms, rmsTotal);
            CONSOLE << iloop << std::endl;

            // Global residual statistics
            Residual resGlobal(res);
            resGlobal.Allgather();
            resGlobal.GetRMS(rms, rmsTotal);
            Residual::BlockIndices maxLocs = resGlobal.MaxLocs();

            //if (iloop == 0 && isub == 0)
            if (iloop < iterStartSER && isub == 0)
            {
                res0 = resGlobal;
                res0.GetRMS(rms0, rmsTotal0);
            }

            if (COMM->MyRank() == 0)
            {
                for (int l = 0; l < 5; ++l)
                {
                    std::cout << std::setw(14) << std::scientific << std::setprecision(5);
                    std::cout << rms[l] / rms0[l];
                }
                std::cout << std::setw(14) << std::scientific << std::setprecision(5);
                std::cout << rmsTotal / rmsTotal0 << std::endl;
                for (int l = 0; l < 6; ++l)
                {
                    std::cout << '\t' << maxLocs[l].BlockID << ":" << maxLocs[l].Index;
                }
                std::cout << std::endl;
            }

            frms << iloop << '\t' << rmsTotal << '\t' << cfl << std::endl;

            bool negative = false;
            std::map<int, std::vector<IndexIJK> > negativeBlockIndices;
            for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
            {
                std::vector<IndexIJK> negativeIndices;
                negative = negative || blocks[iblock]->CheckNegatives(negativeIndices);
                negativeBlockIndices[blocks[iblock]->ID()] = negativeIndices;
            }
            negative = COMM->Any(negative);
            if (negative)
            {
                CONSOLE << "negative quantities detected somewhere, bailing out." << std::endl;
                bailout = true;
                WriteNegativeCells("negatives.vtk", negativeBlockIndices, meshScale);
                //Probe(blocks);
                break;
            }

            if (iloop >= iterStartSER)
                cfl = std::max(cfl0, std::min(cflmax, cfl0 * std::pow(rmsTotal0 / rmsTotal, ser)));
            else
                cfl = cfl0;
            tsEval.SetCFL(cfl);

            // ******************************************
            // Turbulence Solve
            // ******************************************
            if (iterT0 >= 0 && iloop >= iterT0)
            {
                tsEvalT.Evaluate(integratorsT);

                int nistepsT = integratorsT[0]->IntegrationSteps();
                for (int ii = 0; ii < nistepsT; ++ii)
                {
                    // Pre-Integrate Phase
                    for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                    {
                        integratorsT[iblock]->PreIntegrateStart(ii, iteration);
                    }
                    for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                    {
                        integratorsT[iblock]->PreIntegrateFinish(ii, iteration);
                    }

                    // Integrate Phase
                    for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                    {
                        Integrator* integratorT = integratorsT[iblock];
                        integratorT->Integrate(ii, iteration, blocks[iblock]->UT());
                    }
                    if (ii == nistepsT - 1)
                    {
                        for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                        {
                            Block* block = blocks[iblock];
                            IndexRange cr = block->CellRange();
                            Integrator* integratorT = integratorsT[iblock];
                            for (int k = cr.Start.K; k <= cr.End.K; ++k)
                            {
                                for (int j = cr.Start.J; j <= cr.End.J; ++j)
                                {
                                    for (int i = cr.Start.I; i <= cr.End.I; ++i)
                                    {
                                        double* ut = block->UT()(i, j, k);
                                        double* dut = integratorT->DU()(i, j, k);
                                        const double EPS_K = 1.0e-20;
                                        if (ut[0] + dut[0] <= 0.0)
                                            ; //ut[0] = EPS_K;
                                        else
                                            ut[0] += dut[0];
                                        if (ut[1] + dut[1] <= 0.0)
                                            ;
                                        else
                                            ut[1] += dut[1];
                                    }
                                }
                            }
                        }
                    }

                    // Post-Integrate Phase
                    for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                    {
                        integratorsT[iblock]->PostIntegrateStart(ii, iteration);
                    }
                    for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                    {
                        integratorsT[iblock]->PostIntegrateFinish(ii, iteration);
                    }
                    COMM->Barrier();
                }

                // Interface data exchange
                for (InterfaceDataExchangers::iterator iidx = IDXsUT.begin(); iidx != IDXsUT.end(); ++iidx)
                {
                    (*iidx)->Exchange();
                }

                // Compute turbulent eddy viscosity
                for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
                {
                    blocks[iblock]->ComputeTurbulentTransportProperties();
                }

                // Accumulate residuals
                Residual resT = integratorsT[0]->LatestResidual();
                for (size_t iblock = 1; iblock < blocks.size(); ++iblock)
                {
                    resT += integratorsT[iblock]->LatestResidual();
                }
                CONSOLE << "resT:" << resT << std::endl;

                // Global residual
                double rmsT[2], rmsTotalT;
                Residual resGlobalT(resT);
                resGlobalT.Allgather();
                resGlobalT.GetRMS(rmsT, rmsTotalT);
                Residual::BlockIndices maxLocsT = resGlobalT.MaxLocs();

                if (iloop == iterT0)
                {
                    res0T = resGlobalT;
                    res0T.GetRMS(rms0T, rmsTotal0T);
                }

                if (COMM->MyRank() == 0)
                {
                    for (int l = 0; l < 2; ++l)
                    {
                        std::cout << std::setw(14) << std::scientific << std::setprecision(5);
                        std::cout << rmsT[l];
                    }
                    std::cout << std::setw(14) << std::scientific << std::setprecision(5);
                    std::cout << rmsTotalT / rmsTotal0T << std::endl;
                    for (int l = 0; l < 3; ++l)
                    {
                        std::cout << '\t' << maxLocsT[l].BlockID << ":" << maxLocsT[l].Index;
                    }
                    std::cout << std::endl;
                }
            } // end of turbulence solve

            Probe(CONSOLE, blocks, probeRanges);
        } // end of subiteration loop

        if (iloop > 0 && iloop % 100 == 0)
        {
            //WriteSolution("tmp", "tmp_residual", 0, *model, blocks, blockRankMap, integrators);
            std::ostringstream oss;
            oss << "tmp." << std::setw(5) << std::setfill('0') << iloop;
            WriteSolutionCGNS(oss.str().c_str(), blocks, blockRankMap, cgnsSt);
        }

        if (unsteady && iloop % 10 == 0)
        {
            //std::ostringstream oss;
            //oss << "out";
            //WriteSolution(oss.str().c_str(), NULL, iteration.TimeStep(), *model, blocks, blockRankMap, integrators);
            std::ostringstream oss;
            oss << "out." << std::setw(5) << std::setfill('0') << iteration.TimeStep();
            WriteSolutionCGNS(oss.str().c_str(), blocks, blockRankMap, cgnsSt);

            for (size_t iblock = 0; iblock < blocks.size(); ++iblock)
            {
                blocks[iblock]->ShiftTime();
            }
        }
        iteration.AdvanceTimeStep();

        if (bailout)
        {
            break;
        }
    }
    std::cout << "DONE" << std::endl;

    WriteSolutionCGNS("out", blocks, blockRankMap, cgnsSt);
    //WriteSolution("out", "residual", blocks, blockRankMap, integrators);
    //WriteSolution("out", NULL, 0, *model, blocks, blockRankMap, integrators);

    COMM->Finalize();

    return EXIT_SUCCESS;
}

