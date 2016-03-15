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
// $Id: CGNSMeshReader.h 76 2010-10-06 16:50:49Z kato $
#ifndef INCLUDED_CGNS_MESH_READER_H__
#define INCLUDED_CGNS_MESH_READER_H__

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "DomainInfo.h"
#include "Structured.h"

class Block;
class BC;

class CGNSMeshReader
{
public:
    CGNSMeshReader(const char* filename);
    virtual ~CGNSMeshReader();

    const DomainInfo& Domain() const { return mDomain; }
    void ReadMesh(Structured<double>& XYZ, int Z) const;
    void ReadZone(Block** block, int Z);
    void ReadBC(BC** bc, int Z, int BC);
    void ReadPartialMesh(Structured<double>& XYZ, int Z, const IndexRange& range) const;

protected:
    void ReadBlockInfo();
    void GetZone(int& fn, int B, int Z) const;
    void Close(int fn) const;

private:
    std::string mFileName;
    DomainInfo mDomain;
};

#endif // INCLUDED_CGNS_MESH_READER_H__

