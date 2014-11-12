// $Id: RampMesh.h 4 2010-03-06 14:10:00Z kato $
#ifndef INCLUDED_RAMP_MESH_H__
#define INCLUDED_RAMP_MESH_H__

#include "Structured.h"

class RampMesh
{
public:
    RampMesh() {}

    void GenerateMesh(Structured<double>& xyz, double LX, double LY, double LZ, double LX1, double delta1) const;

protected:

private:
};

#endif // INCLUDED_RAMP_MESH_H__

