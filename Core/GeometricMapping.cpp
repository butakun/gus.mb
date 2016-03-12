#include "GeometricMapping.h"
#include "ANN/ANN.h"

GeometricMapping::GeometricMapping(const BlockPatches& patches1, const BlockPatches& patches2)
:   mPatches1(patches1), mPatches2(patches2)
{
}

GeometricMapping::~GeometricMapping()
{
}

void
GeometricMapping::InitializeMapping()
{
}

void
GeometricMapping::GenerateMapping(const Vector3& angle1, const Vector3& angle2)
{
}
