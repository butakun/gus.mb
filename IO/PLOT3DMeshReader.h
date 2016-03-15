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
// $Id: PLOT3DMeshReader.h 4 2010-03-06 14:10:00Z kato $
#ifndef INCLDUED_HIRO_PLOT3D_MESH_READER_H__
#define INCLDUED_HIRO_PLOT3D_MESH_READER_H__

#include "MeshReader.h"
#include <string>

class PLOT3DMeshReader : public MeshReader
{
public:
    PLOT3DMeshReader(const char* filename);
    virtual ~PLOT3DMeshReader();

    virtual void Read(Structured<double>& xyz, IndexRange& range);

protected:

private:
    std::string mFileName;
};

#endif // INCLDUED_HIRO_PLOT3D_MESH_READER_H__

