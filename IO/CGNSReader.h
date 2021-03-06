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
// $Id: CGNSReader.h 296 2013-08-30 06:53:12Z kato $
#ifndef INCLUDED_CGNS_READER_H__
#define INCLUDED_CGNS_READER_H__

class Block;
#include "Structured.h"
#include "Vector3.h"
#include "Physics.h"
#include "RigidBodyMotion.h"
#include <string>
#include <algorithm>

class CGNSStructure
{
public:
    class Zone
    {
    public:
        Zone(int z, const char* name, const IndexIJK& vsize, const IndexIJK& csize)
        : Z(z), Name(name), VertexSize(vsize), CellSize(csize), Motion(NULL) {}
        int Z;
        std::string Name;
        IndexIJK VertexSize; // # of vertices in 3 directions
        IndexIJK CellSize; // # of cells in 3 directions
        RigidBodyMotion* Motion;

        static bool equal(const Zone& zone, int z) { return zone.Z == z; }
        bool equal_to(int z) const { return this->Z == z; }
    };
    typedef std::vector<Zone> Zones_t;

    class Base
    {
    public:
        Base(int b, const char* name, int cdim, int pdim)
        : B(b), Name(name), CellDim(cdim), PhysDim(pdim) {}
        int B;
        std::string Name;
        int CellDim, PhysDim;
        Zones_t Zones;

        Zone& FindZone(int Z)
        {
            Zones_t::iterator i = std::find_if(Zones.begin(), Zones.end(), std::bind2nd(std::mem_fun_ref(&Zone::equal_to), Z));
            return *i;
        }
    };
    typedef std::vector<Base> Bases_t;

    Bases_t Bases;

    CGNSStructure() {}

protected:

private:
};

class CGNSReader
{
public:
    CGNSReader(const char* filename);
    virtual ~CGNSReader();

    CGNSStructure ReadStructure() const;

    void ReadMesh(int zone, Structured<double>& XYZ, std::string& zoneName, const IndexRange rangeToRead = IndexRange(0, 0, 0, 0, 0, 0));
    void ReadFlowSolution(int zone, const Block& block, Structured<double>& U, const Physics& phys);
    void ReadTurbulenceSolution(int zone, const Block& block, const Structured<double>& U, Structured<double>& UT, const Physics& phy, const std::string& model);

    //void ReadFlowSolutionDim(int zone, Structur

protected:

private:
    int fn;
};

#endif // INCLUDED_CGNS_READER_H__

