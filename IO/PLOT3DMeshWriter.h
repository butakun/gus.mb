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
// $Id: PLOT3DMeshWriter.h 76 2010-10-06 16:50:49Z kato $
#ifndef INCLUDED_PLOT3D_MESH_WRITER_H__
#define INCLUDED_PLOT3D_MESH_WRITER_H__

#include <vector>

class Block;

class PLOT3DMeshWriter
{
public:
    PLOT3DMeshWriter() {}

    void AddBlock(const Block& block);
    void Write(const char* filename) const;

    void Write(const char* filename, const Block& block) const;

protected:

private:
    std::vector<const Block*> mBlocks;
};

#endif // INCLUDED_PLOT3D_MESH_WRITER_H__

