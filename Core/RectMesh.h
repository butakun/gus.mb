// $Id: RectMesh.h 55 2010-07-27 14:21:56Z kato $
#ifndef INCLUDED_RECT_MESH_H__
#define INCLUDED_RECT_MESH_H__

#include "Structured.h"

class RectMesh
{
public:
    RectMesh() {}

    /// expects xyz is already allocated
    void GenerateMesh(Structured<double>& xyz, double LX, double LY, double LZ, double dywall = -1.0) const;

protected:

private:
};

#endif // INCLUDED_RECT_MESH_H__

