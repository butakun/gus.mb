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
#include "Physics.h"
#include "Roster.h"
#include "FlowModel.h"
#include "Block.h"
#include "CylinderMesh.h"
#include "Connectivity1to1.h"
#include "StructuredDataExchanger.h"
#include "VTKWriter.h"
#include <cstdlib>
#include <iostream>

void Dump(
    const Structured<double>& U1,
    const Structured<double>& U2,
    const IndexRange& cellRange
    )
{
    int k = 1;
    for (int j = cellRange.Start.J; j <= cellRange.End.J; ++j)
    {
        for (int i = cellRange.Start.I; i <= cellRange.End.I; ++i)
        {
            std::cout << int(U1(i, j, k)[0]) << ' ';
        }
        std::cout << '|';
        for (int i = cellRange.Start.I; i <= cellRange.End.I; ++i)
        {
            std::cout << int(U2(i, j, k)[0]) << ' ';
        }
        std::cout << std::endl;
    }
}

void Step(Structured<double>& U, const IndexRange& cr)
{
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.End.I; i >= cr.Start.I; --i)
            {
                for (int l = 0; l < U.DOF(); ++l)
                {
                    U(i, j, k)[l] = U(i - 1, j, k)[l];
                }
            }
        }
    }
}

int main(int argc, char** argv)
{
    Communicator::Initialize(&argc, &argv);

    FlowModel model;

    int imax = 1, jmax = 1, kmax = 1;
    IndexRange meshRange(0, 0, 0, imax, jmax, kmax);
    std::cout << meshRange.Size() << std::endl;

    Block* block1 = Block::New(1, meshRange);
    Block* block2 = Block::New(2, meshRange);
    CylinderMesh mesher;
    double dTheta = 360.0 / 15;
    mesher.GenerateMesh(block1->XYZ(), 0.0, 1.0, 1.0, 3.0, 0.0, dTheta);
    mesher.GenerateMesh(block2->XYZ(), 0.0, 1.0, 1.0, 3.0, dTheta, dTheta);
    IndexRange cr1 = block1->CellRange();
    IndexRange cr2 = block2->CellRange();

    int tag = 1234;
    int t[3] = {1, 2, 3};
    Connectivity1to1 conn1(
        IndexRange(0, jmax, 0, imax, jmax, kmax), JNEG, *block1,
        IndexTransform(t, IndexIJK(0, jmax, 0), IndexIJK(0, 0, 0)),
        block2->ID(), tag
        );
    block1->RegisterBC(&conn1);
    std::cout << "debug" << std::endl;

    Connectivity1to1 conn2(
        IndexRange(0, 0, 0, imax, 0, kmax), J, *block2,
        IndexTransform(t, IndexIJK(0, 0, 0), IndexIJK(0, jmax, 0)),
        block1->ID(), tag
        );
    block2->RegisterBC(&conn2);

    for (int ib = 1; ib <= 2; ++ib)
    {
        Block* block = dynamic_cast<Block*>(Roster::GetInstance()->GetBlock(ib));
        assert(block != NULL);
        Structured<double>& XYZ = block->XYZ();
        Structured<double>& U = block->U();
        IndexRange cr = block->CellRange();
        IndexIterator itor(cr);
        while (!itor.IsEnd())
        {
            IndexIJK ijk = itor.Index();
            Vector3 p1(XYZ(ijk + IndexIJK(-1, -1, -1)));
            Vector3 p2(XYZ(ijk + IndexIJK( 0, -1, -1)));
            Vector3 p3(XYZ(ijk + IndexIJK( 0,  0, -1)));
            Vector3 p4(XYZ(ijk + IndexIJK(-1,  0, -1)));
            Vector3 p5(XYZ(ijk + IndexIJK(-1, -1,  0)));
            Vector3 p6(XYZ(ijk + IndexIJK( 0, -1,  0)));
            Vector3 p7(XYZ(ijk + IndexIJK( 0,  0,  0)));
            Vector3 p8(XYZ(ijk + IndexIJK(-1,  0,  0)));
            Vector3 pc = 0.125 * (p1 + p2 + p3 + p4);
            Vector3 vr(0.0, pc.Y(), pc.Z());
            vr.Normalize();
            U(ijk)[0] = 1.0;
            U(ijk)[1] = 1.0 * vr.X();
            U(ijk)[2] = 1.0 * vr.Y();
            U(ijk)[3] = 1.0 * vr.Z();
            U(ijk)[4] = 1.0;
            std::cout << "Block " << block->ID() << " " << ijk << std::endl;
            itor.Advance();
        }
    }

    StructuredDataExchanger ex1(*block1, &model, block1->U());
    StructuredDataExchanger ex2(*block2, &model, block2->U());

    ex1.Start();
    ex2.Start();
    ex1.Finish();
    ex2.Finish();

#if 0
    Dump(block1->U(), block2->U(), cr1);

    for (int step = 1; step <= 10; ++step)
    {
        Step(block1->U(), cr1);
        Step(block2->U(), cr2);
        ex1.Start();
        ex2.Start();
        ex1.Finish();
        ex2.Finish();

        Dump(block1->U(), block2->U(), cr1);
    }
#endif

    VTKWriter* writer;
    writer = new VTKWriter("test_perio.1.vts", VTKWriter::XML);
    writer->AddMesh(block1->XYZ());
    writer->Write();
    delete writer;
    writer = new VTKWriter("test_perio.2.vts", VTKWriter::XML);
    writer->AddMesh(block2->XYZ());
    writer->Write();
    delete writer;
    writer = new VTKWriter("dummy", VTKWriter::XML);
    std::vector<std::string> filenames;
    filenames.push_back("test_perio.1.vts");
    filenames.push_back("test_perio.2.vts");
    writer->Finalize("test_perio.vtm", filenames);
    delete writer;

    return EXIT_SUCCESS;
}

