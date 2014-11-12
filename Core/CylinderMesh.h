// $Id: CylinderMesh.h 41 2010-07-06 10:48:07Z kato $
#ifndef INCLUDED_CYLINDER_MESH_H__
#define INCLUDED_CYLINDER_MESH_H__

#include "Structured.h"

class CylinderMesh
{
public:
    CylinderMesh() {}

    void GenerateMesh(Structured<double>& xyz, double X0, double LX, double RHub, double RShroud, double theta0, double dTheta) const;

protected:

private:
};

#endif // INCLUDED_CYLINDER_MESH_H__

