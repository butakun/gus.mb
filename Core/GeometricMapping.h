#ifndef INCLUDED_GEOMETRIC_MAPPING_H__
#define INCLUDED_GEOMETRIC_MAPPING_H__

#include "MeshMapper.h"
#include "BlockPatch.h"
#include "Vector3.h"

class GeometricMapping : public MeshMapper
{
public:
    GeometricMapping(const BlockPatches& patches1, const BlockPatches& patches2);
    virtual ~GeometricMapping();

    virtual void InitializeMapping();
    virtual void GenerateMapping(const Vector3& angle1, const Vector3& angle2);

protected:

private:
    BlockPatches mPatches1, mPatches2;
};

#endif // INCLUDED_GEOMETRIC_MAPPING_H__

