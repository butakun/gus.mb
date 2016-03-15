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
#ifndef INCLUDED_GEOMETRY_TOOLS_H__
#define INCLUDED_GEOMETRY_TOOLS_H__

#include "Vector3.h"

namespace GeometryTools
{

// approximate area of a 3-d quadrilateral, sum of two triangles (123 and 134)
// 4 --- 3
// |     |
// 1 --- 2
inline
double QuadArea(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4)
{
    return 0.5 * (cross_product(p2 - p1, p3 - p1).Mag() + cross_product(p3 - p1, p4 - p1).Mag());
}

}

#endif // INCLUDED_GEOMETRY_TOOLS_H__

