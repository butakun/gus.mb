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
// $Id: TestCGNSMeshReader.cpp 19 2010-05-14 15:28:18Z kato $

#include "Block.h"
#include "CGNSMeshReader.h"
#include "PLOT3DMeshWriter.h"
#include <iostream>

int
main(int argc, char** argv)
{
    Block* block1;
    Block* block2;
    CGNSMeshReader reader(argv[1]);

    for (int Z = 1; Z <= reader.Domain().GetBlockInfos().size(); ++Z)
    {
        const BlockInfo& bi = reader.Domain().FindBlockInfo(Z);
        std::cout << "Zone " << Z << std::endl;
        std::cout << bi << std::endl;
    }

#if 0
    reader.ReadZone(&block1, 1);
    reader.ReadZone(&block2, 2);
#endif

#if 0
    std::cout << "XYZ" << std::endl;
    for (int k = 0; k < 10; ++k)
    {
        for (int j = 0; j < 10; ++j)
        {
            for (int i = 0; i < 10; ++i)
            {
                double* xyz = block1->XYZ()(i, j, k);
                std::cout << xyz[0] << ' ' << xyz[1] << ' ' << xyz[2] << std::endl;
            }
        }
    }
    return EXIT_SUCCESS;
#endif

#if 0
    PLOT3DMeshWriter writer;
    writer.AddBlock(*block1);
    writer.AddBlock(*block2);
    writer.Write("out.xyz");
#endif

    return EXIT_SUCCESS;
}

