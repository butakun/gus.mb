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
// $Id: Wall.h 37 2010-06-11 16:00:18Z kato $
#ifndef INCLUDED_WALL_H__
#define INCLUDED_WALL_H__

#include "Structured.h"

class Wall
{
public:
    Wall(const IndexRange& meshRange)
    :   mXYZ(3, meshRange)
    {}

    ~Wall()
    {
        delete[] mXYZ.Data;
    }

    const Structured<double>& XYZ() const { return mXYZ; }

protected:

private:
    Structured<double> mXYZ;
};

#endif // INCLUDED_WALL_H__

