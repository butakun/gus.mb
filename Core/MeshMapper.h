// $Id: MeshMapper.h 261 2013-01-05 12:01:20Z kato $
#ifndef INCLUDED_MESH_MAPPER_H__
#define INCLUDED_MESH_MAPPER_H__

#include "Structured.h"
#include "Model.h"

class Vector3;
class BlockPatch;

class MeshMapper
{
public:
    virtual ~MeshMapper() {}

    virtual void InitializeMapping() = 0;
    virtual void MapOn(const BlockPatch& bpd, const Structured<double>& XYZDonor, const Vector3& angleSelf, const Vector3& angleDonor, bool debug = false) = 0;
    virtual bool AllMapped() const = 0;

    virtual void MapData(const Model& model, Structured<double>& U, const BlockPatch& bpd, const Structured<double>& UDonor) const = 0;

protected:

private:
};

#endif // INCLUDED_MESH_MAPPER_H__

