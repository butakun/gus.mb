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

