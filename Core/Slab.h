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
// $Id: Slab.h 41 2010-07-06 10:48:07Z kato $
#ifndef INCLUDED_SLAB_H__
#define INCLUDED_SLAB_H__

#include "Structured.h"
#include <string>
#include <map>

class Slab
{
public:
    Slab(int blockID, const char* name);
    virtual ~Slab();

    int BlockID() const { return mBlockID; }
    const std::string& Name() const { return mName; }

    Structured<double>& AddData(const char* name, int dof, const IndexRange& range);
    Structured<double>& GetData(const char* name);

protected:

private:
    typedef std::map<std::string, Structured<double> > DataMap;

    int mBlockID;
    std::string mName;
    DataMap mDataMap;
};

#endif // INCLUDED_SLAB_H__

