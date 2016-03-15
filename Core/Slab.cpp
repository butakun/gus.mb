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
// $Id: Slab.cpp 41 2010-07-06 10:48:07Z kato $

#include "Slab.h"
#include <cassert>

Slab::Slab(int blockID, const char* name)
:   mBlockID(blockID), mName(name)
{
}

Slab::~Slab()
{
    for (DataMap::iterator i = mDataMap.begin(); i != mDataMap.end(); ++i)
    {
        Structured<double>& d = i->second;
        delete[] d.Data;
    }
}

Structured<double>&
Slab::AddData(const char* name, int dof, const IndexRange& range)
{
    DataMap::const_iterator i = mDataMap.find(name);
    assert(i == mDataMap.end());

    Structured<double> d(dof, range);
    mDataMap[name] = d;

    return mDataMap[name];
}

Structured<double>&
Slab::GetData(const char* name)
{
    DataMap::iterator i = mDataMap.find(name);
    assert(i != mDataMap.end());

    return i->second;
}

