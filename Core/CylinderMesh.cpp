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
// $Id: CylinderMesh.cpp 41 2010-07-06 10:48:07Z kato $

#include "CylinderMesh.h"
#include "Vector3.h"
#include <cmath>
#include <cassert>

void
CylinderMesh::GenerateMesh(Structured<double>& xyz, double X0, double LX, double RHub, double RShroud, double theta0, double dTheta) const
{
    assert(xyz.Data != NULL);

    IndexRange range = xyz.GetRange();

    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        double zeta = double(k - range.Start.K) / double(range.End.K - range.Start.K);
        double R = (1.0 - zeta) * RHub + zeta * RShroud;
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            double eta = double(j - range.Start.J) / double(range.End.J - range.Start.J);
            double theta = (1.0 - eta) * theta0 + eta * (theta0 + dTheta);
            double y, z;
            y = R * std::sin(theta * M_PI / 180.0);
            z = R * std::cos(theta * M_PI / 180.0);
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                double xi = double(i - range.Start.I) / double(range.End.I - range.Start.I);
                double x = (1.0 - xi) * X0 + xi * (X0 + LX);
                xyz(i, j, k)[0] = x;
                xyz(i, j, k)[1] = y;
                xyz(i, j, k)[2] = z;
            }
        }
    }
}

