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
#include "Block.h"
#include "RectMesh.h"
#include "Connectivity1to1.h"
#include "StructuredDataExchanger.h"
#include "FlowModel.h"
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

    int imax = 10, jmax = 10, kmax = 1;
    IndexRange meshRange(0, 0, 0, imax, jmax, kmax);
    std::cout << meshRange.Size() << std::endl;

    Block* block1 = Block::New(1, meshRange);
    Block* block2 = Block::New(2, meshRange);
    RectMesh mesher;
    mesher.GenerateMesh(block1->XYZ(), 1.0, 1.0, 0.1);
    mesher.GenerateMesh(block2->XYZ(), 1.0, 1.0, 0.1);
    IndexRange cr1 = block1->CellRange();
    IndexRange cr2 = block2->CellRange();

    for (int k = meshRange.Start.K; k <= meshRange.End.K; ++k)
    {
        for (int j = meshRange.Start.J; j <= meshRange.End.J; ++j)
        {
            for (int i = meshRange.Start.I; i <= meshRange.End.I; ++i)
            {
                block2->XYZ()(i, j, k)[0] += 1.0;
            }
        }
    }

    int tag = 1234;
    int t[3] = {1, 2, 3};
    Connectivity1to1 conn1(
        IndexRange(imax, 0, 0, imax, jmax, kmax), INEG, *block1,
        IndexTransform(t, IndexIJK(imax, 0, 0), IndexIJK(0, 0, 0)),
        block2->ID(), tag
        );
    block1->RegisterBC(&conn1);
    std::cout << "debug" << std::endl;

    Connectivity1to1 conn2(
        IndexRange(0, 0, 0, 0, jmax, kmax), I, *block2,
        IndexTransform(t, IndexIJK(0, 0, 0), IndexIJK(imax, 0, 0)),
        block1->ID(), tag
        );
    block2->RegisterBC(&conn2);

    double U0[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    block1->U().SetTo(U0);
    block2->U().SetTo(U0);
    for (int i = 1; i <= imax; ++i)
    {
        block1->U()(i, jmax / 2, 1)[0] = 1.0;
    }

    StructuredDataExchanger ex1(*block1, &model, block1->U());
    StructuredDataExchanger ex2(*block2, &model, block2->U());

    ex1.Start();
    ex2.Start();
    ex1.Finish();
    ex2.Finish();

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

    return EXIT_SUCCESS;
}

