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
// $Id: BCEnforce2D.cpp 175 2012-01-04 06:01:45Z kato $

#include "BCEnforce2D.h"

BCEnforce2D::BCEnforce2D(Type type)
:   BC(IndexRange(-1, -1, -1, -1, -1, -1)), mType(type)
{
}

BCEnforce2D::~BCEnforce2D()
{
}

void
BCEnforce2D::Apply(const Block& block, Structured<double>& U)
{
    IndexRange r = U.GetRange();
    for (IndexIterator it(r); !it.IsEnd(); it.Advance())
    {
        IndexIJK ijk = it.Index();

        // FIXME: only TWOD_Z is implemented.
        double* UU = U(ijk);
        Vector3 e(UU[1], UU[2], 0.0);
        e.Normalize();
        double rhoV = std::sqrt(UU[1] * UU[1] + UU[2] * UU[2] + UU[3] * UU[3]);
        UU[1] = rhoV * e.X();
        UU[2] = rhoV * e.Y();
        UU[3] = rhoV * e.Z();
    }
}

void
BCEnforce2D::ApplyTurb(const Block& block, Structured<double>& UT)
{
}

