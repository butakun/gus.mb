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

