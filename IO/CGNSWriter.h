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
// $Id: CGNSWriter.h 294 2013-08-23 14:11:39Z kato $
#ifndef INCLUDE_CGNS_WRITER_H__
#define INCLUDE_CGNS_WRITER_H__

#include "Structured.h"
#include "Block.h"
#include "Physics.h"
#include <string>

class CGNSStructure;

class CGNSWriter
{
public:
    CGNSWriter(const char* filename, bool modify = true);
    ~CGNSWriter();

    void WriteStructure(const CGNSStructure& s);
    void WriteFlowSolution(int zone, const Block& block, const Structured<double>& U, const Physics& phys);
    void WriteTurbulenceSolution(int zone, const Block& block, const Structured<double>& U, const Structured<double>& UT, const Physics& phys, const std::string& model);

    int WriteBase();
    int WriteZone(int B, const IndexRange& meshRange, const char* name = NULL);

protected:

private:
    int fn;
};

#endif // INCLUDE_CGNS_WRITER_H__

