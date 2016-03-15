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
// $Id: PLOT3DMeshWriter.cpp 89 2011-01-28 14:16:26Z kato $

#include "PLOT3DMeshWriter.h"
#include "Block.h"
#include <fstream>

void
PLOT3DMeshWriter::AddBlock(const Block& block)
{
    mBlocks.push_back(&block);
}

void WriteBlock(std::ostream& o, const Block& block)
{
    IndexRange mr = block.MeshRange();
    const Structured<double>& xyz = block.XYZ();

    for (int l = 0; l < 3; ++l)
    {
        for (int k = mr.Start.K; k <= mr.End.K; ++k)
        {
            for (int j = mr.Start.J; j <= mr.End.J; ++j)
            {
                for (int i = mr.Start.I; i <= mr.End.I; ++i)
                {
                    o << xyz(i, j, k)[l] << std::endl;
                }
            }
        }
    }
}

void
PLOT3DMeshWriter::Write(const char* filename) const
{
    std::ofstream f(filename);

    f << mBlocks.size() << std::endl;
    for (size_t i = 0; i < mBlocks.size(); ++i)
    {
        const Block& block = *mBlocks[i];
        IndexIJK shape = block.MeshRange().Shape();
        f << shape.I << '\t' << shape.J << '\t' << shape.K << std::endl;
    }
    for (size_t i = 0; i < mBlocks.size(); ++i)
    {
        const Block& block = *mBlocks[i];
        WriteBlock(f, block);
    }

    f.close();
}

void
PLOT3DMeshWriter::Write(const char* filename, const Block& block) const
{
    std::ofstream f(filename);

    IndexIJK shape = block.MeshRange().Shape();
    f << shape.I << '\t' << shape.J << '\t' << shape.K << std::endl;
    WriteBlock(f, block);

    f.close();
}

