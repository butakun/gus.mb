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
// $Id: Model.h 253 2012-07-20 10:03:41Z kato $
#ifndef INCLUDED_MODEL_H__
#define INCLUDED_MODEL_H__

#include "Block.h"

class Model
{
public:
    virtual ~Model() {}

    // Coordinate frame transform
    virtual void FromGlobalToLocal(double* ULocal, double* UGlobal, const Block& block, const IndexIJK& ijk) const = 0;
    virtual void FromLocalToGlobal(double* UGlobal, double* ULocal, const Block& block, const IndexIJK& ijk) const = 0;
};

#endif // INCLUDED_MODEL_H__

