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

