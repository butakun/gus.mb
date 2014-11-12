// $Id: BLMesh.h 121 2011-08-18 16:08:26Z kato $
#ifndef INCLUDED_BL_MESH_H__
#define INCLUDED_BL_MESH_H__

#include "Structured.h"

class BLMesh
{
public:
    BLMesh() {}

    /// expects xyz is already allocated
    void GenerateMesh(Structured<double>& xyz, double LX, double LY, double LZ, int NLE, double XLE, double dywall = -1.0) const;

protected:

private:
};

#endif // INCLUDED_BL_MESH_H__

