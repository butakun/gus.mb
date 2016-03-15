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
// $Id: ConvertCGNSToPLOT3D.cpp 206 2012-03-07 03:42:14Z kato $

#include "CGNSReader.h"
#include "PLOT3DMeshWriter.h"
#include "Structured.h"
#include "Block.h"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstdlib>

void ConvertZone(const char* cgnsFileName, const char* plot3dFileName, int Z)
{
    CGNSReader reader(cgnsFileName);

    Structured<double> XYZ;
    std::string zoneName;
    reader.ReadMesh(Z, XYZ, zoneName);

    Block* block = Block::New(Z, XYZ.GetRange());
    block->XYZ().SetTo(XYZ);

    std::ostringstream oss;
    oss << plot3dFileName << "." << Z << ".xyz";
    PLOT3DMeshWriter meshWriter;
    meshWriter.Write(oss.str().c_str(), *block);
}

int main(int argc, char** argv)
{
    std::cout << argc << std::endl;
    if (argc < 3)
    {
        std::cerr << argv[0] << " CGNSFileName Zone1 Zone2 Zone3 ..." << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<int> zones;
    for (int i = 2; i < argc; ++i)
    {
        std::istringstream iss(argv[i]);
        int z;
        iss >> z;
        zones.push_back(z);
    }

    const char* cgnsFileName = argv[1];
    std::string plot3dFileName = "plot3d";

    for (std::vector<int>::const_iterator i = zones.begin(); i != zones.end(); ++i)
    {
        int Z = *i;
        ConvertZone(cgnsFileName, plot3dFileName.c_str(), Z);
    }

    return EXIT_SUCCESS;
}

