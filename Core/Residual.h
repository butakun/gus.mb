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

