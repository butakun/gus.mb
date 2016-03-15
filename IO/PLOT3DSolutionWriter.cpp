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
// $Id: PLOT3DSolutionWriter.cpp 79 2010-12-09 13:51:48Z kato $

#include "PLOT3DSolutionWriter.h"
#include <fstream>

PLOT3DSolutionWriter::PLOT3DSolutionWriter(const char* name)
:   mFileName(name)
{
}

void
PLOT3DSolutionWriter::Write(const Structured<double>& U, const IndexRange& meshRange)
{
    double FSMach = 1.0, Alpha = 0.0, Rey = 1.0e5;

    std::ofstream f(mFileName.c_str());

    IndexIJK shape = meshRange.Shape();

    f << shape.I << '\t' << shape.J << '\t' << shape.K << std::endl;
    f << FSMach << '\t' << Alpha << '\t' << Rey << '\t' << 0.0 << std::endl;

    for (int l = 0; l < 5; ++l)
    {
        for (int k = meshRange.Start.K; k <= meshRange.End.K; ++k)
        {
            for (int j = meshRange.Start.J; j <= meshRange.End.J; ++j)
            {
                for (int i = meshRange.Start.I; i <= meshRange.End.I; ++i)
                {
                    IndexIJK i1(i, j, k), i2(i + 1, j, k), i3(i + 1, j + 1, k), i4(i, j + 1, k),
                        i5(i, j, k + 1), i6(i + 1, j, k + 1), i7(i + 1, j + 1, k + 1), i8(i, j + 1, k + 1);
                    double u;
                    u = 0.125 * (U(i1)[l] + U(i2)[l] + U(i3)[l] + U(i4)[l] + U(i5)[l] + U(i6)[l] + U(i7)[l] + U(i8)[l]);
                    f << u << std::endl;
                }
            }
        }
    }
}

