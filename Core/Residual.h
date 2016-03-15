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
// $Id: Residual.h 154 2011-10-14 05:11:44Z kato $
#ifndef INCLUDED_RESIDUAL_H__
#define INCLUDED_RESIDUAL_H__

#include "Structured.h"
#include <iostream>
#include <vector>

class IndexIJK;

class Residual
{
public:
    class BlockIndex;
    typedef std::vector<BlockIndex> BlockIndices;

    Residual();
    Residual(int dof, int numPoints);
    Residual(int dof, double* square, double* maxSquare, const BlockIndices& maxlocs, int numPoints);
    Residual(const Structured<double>& R, const IndexRange& range, int blockID);
    Residual(const Residual& r);
    ~Residual();

    Residual& operator = (const Residual& r);

    int DOF() const { return mDOF; }
    void GetRMS(double* rms, double& rmsTotal) const;
    const BlockIndices& MaxLocs() const { return mMaxLocs; }
    double* Square() const { return mSquare; }

    void Allgather();

    std::ostream& Dump(std::ostream& o) const;

    Residual operator + (const Residual& r) const;
    Residual& operator += (const Residual& r) { *this = *this + r; return *this; }

    class BlockIndex
    {
    public:
        BlockIndex()
        : BlockID(-1), Index(0, 0, 0)
        {}

        BlockIndex(int blockID, const IndexIJK& index)
        : BlockID(blockID), Index(index)
        {}

        int BlockID;
        IndexIJK Index;
    };

protected:

private:
    int mDOF;
    double* mSquare;        // [DOF]
    double* mMaxSquare;     // [DOF + 1]
    BlockIndices mMaxLocs;  // [DOF + 1]
    int mNumPoints;
};

inline
std::ostream& operator << (std::ostream& o, const Residual& r)
{
    r.Dump(o);
    return o;
}

inline
std::ostream& operator << (std::ostream& o, const Residual::BlockIndex& bi)
{
    o << bi.BlockID << ':' << bi.Index;
    return o;
}

#endif // INCLUDED_RESIDUAL_H__

