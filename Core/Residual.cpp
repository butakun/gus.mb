// $Id: Residual.cpp 180 2012-01-13 08:48:04Z kato $

#include "Communicator.h"
#include "Residual.h"
#include <cstdlib>
#include <mpi.h>

Residual::Residual()
:   mDOF(0), mSquare(NULL), mMaxSquare(NULL), mMaxLocs(), mNumPoints(0)
{
}

Residual::Residual(int dof, int numPoints)
:   mDOF(dof), mSquare(new double[mDOF]), mMaxSquare(new double[mDOF + 1]), mMaxLocs(mDOF + 1),
    mNumPoints(numPoints)
{
    for (int i = 0; i < mDOF; ++i)
    {
        mSquare[i] = 0.0;
        mMaxSquare[i] = 0.0;
        mMaxLocs[i] = BlockIndex(0, IndexIJK(0, 0, 0));
    }
    mMaxSquare[mDOF] = 0.0;
    mMaxLocs[mDOF] = BlockIndex(0, IndexIJK(0, 0, 0));
}

Residual::Residual(int dof, double* square, double* maxSquare, const BlockIndices& maxlocs, int numPoints)
:   mDOF(dof), mSquare(new double[mDOF]), mMaxSquare(new double[mDOF + 1]), mMaxLocs(maxlocs),
    mNumPoints(numPoints)
{
    assert(mMaxLocs.size() == size_t(mDOF + 1));

    for (int i = 0; i < mDOF; ++i)
    {
        mSquare[i] = square[i];
        mMaxSquare[i] = maxSquare[i];
    }
}

Residual::Residual(const Residual& r)
:   mDOF(r.mDOF), mSquare(new double[r.mDOF]), mMaxSquare(new double[r.mDOF + 1]), mMaxLocs(r.mMaxLocs),
    mNumPoints(r.mNumPoints)
{
    assert(mMaxLocs.size() == size_t(mDOF + 1));

    for (int i = 0; i < mDOF; ++i)
    {
        mSquare[i] = r.mSquare[i];
        mMaxSquare[i] = r.mMaxSquare[i];
    }
}

Residual::Residual(const Structured<double>& R, const IndexRange& range, int blockID)
:   mDOF(R.DOF()), mSquare(new double[mDOF]), mMaxSquare(new double[mDOF + 1]), mMaxLocs(mDOF + 1),
    mNumPoints(range.Count())
{
    std::vector<IndexIJK> maxindices;
    R.ReduceSquared(mSquare, mMaxSquare, maxindices, range);
    for (int i = 0; i < mDOF + 1; ++i)
    {
        mMaxLocs[i] = BlockIndex(blockID, maxindices[i]);
    }
}

Residual::~Residual()
{
    delete[] mSquare;
    delete[] mMaxSquare;
}

void
Residual::GetRMS(double* rms, double& rmsTotal) const
{
    rmsTotal = 0.0;
    for (int i = 0; i < mDOF; ++i)
    {
        rms[i] = std::sqrt(mSquare[i] / double(mNumPoints));
        rmsTotal += mSquare[i];
    }
    rmsTotal = std::sqrt(rmsTotal / double(mNumPoints) / double(mDOF));
}

void
Residual::Allgather()
{
    Communicator* COMM = Communicator::GetInstance();
    size_t size = COMM->Size();

    double* allSquares = new double[size * mDOF];
    double* allMaxSquares = new double[size * (mDOF + 1)];
    int* allMaxLocs = new int[size * (mDOF + 1) * 4]; // [BlockID, I, J, K] for each maxloc
    int* allNumPoints = new int[size];

    int* maxlocs = new int[(mDOF + 1) * 4];
    for (int i = 0; i < mDOF + 1; ++i)
    {
        maxlocs[i * 4 + 0] = mMaxLocs[i].BlockID;
        maxlocs[i * 4 + 1] = mMaxLocs[i].Index.I;
        maxlocs[i * 4 + 2] = mMaxLocs[i].Index.J;
        maxlocs[i * 4 + 3] = mMaxLocs[i].Index.K;
    }

    int err;
    err = MPI_Allgather(mSquare, mDOF, MPI_DOUBLE, allSquares, mDOF, MPI_DOUBLE, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
    err = MPI_Allgather(mMaxSquare, mDOF + 1, MPI_DOUBLE, allMaxSquares, mDOF + 1, MPI_DOUBLE, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
    err = MPI_Allgather(maxlocs, (mDOF + 1) * 4, MPI_INT, allMaxLocs, (mDOF + 1) * 4, MPI_INT, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
    err = MPI_Allgather(&mNumPoints, 1, MPI_INT, allNumPoints, 1, MPI_INT, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);

    for (int i = 0; i < mDOF; ++i)
    {
        mSquare[i] = 0.0;
    }
    mNumPoints = 0;

    for (size_t i = 0; i < size; ++i)
    {
        for (int j = 0; j < mDOF; ++j)
        {
            mSquare[j] += allSquares[mDOF * i + j];
        }
        mNumPoints += allNumPoints[i];

        int* maxloc = &allMaxLocs[4 * (mDOF + 1) * i];
        for (int j = 0; j < mDOF + 1; ++j)
        {
            if (mMaxSquare[j] < allMaxSquares[(mDOF + 1) * i + j])
            {
                mMaxSquare[j] = allMaxSquares[(mDOF + 1) * i + j];
                BlockIndex maxlocBlockIndex(maxloc[0], IndexIJK(maxloc[1], maxloc[2], maxloc[3]));
                mMaxLocs[j] = maxlocBlockIndex;
            }
        }
    }

    delete[] allSquares;
    delete[] allMaxSquares;
    delete[] allMaxLocs;
    delete[] allNumPoints;
    delete[] maxlocs;
}

Residual&
Residual::operator = (const Residual& r)
{
    if (mSquare != NULL)
        delete[] mSquare;
    if (mMaxSquare != NULL)
        delete[] mMaxSquare;

    mDOF = r.mDOF;
    mSquare = new double[mDOF];
    mMaxSquare = new double[mDOF + 1];
    mMaxLocs.resize(r.mMaxLocs.size());
    for (int i = 0; i < mDOF; ++i)
    {
        mSquare[i] = r.mSquare[i];
    }
    for (int i = 0; i < mDOF + 1; ++i)
    {
        mMaxSquare[i] = r.mMaxSquare[i];
        mMaxLocs[i] = r.mMaxLocs[i];
    }
    mNumPoints = r.mNumPoints;

    return *this;
}

std::ostream&
Residual::Dump(std::ostream& o) const
{
    double* rms = new double[mDOF];
    double rmsTotal;
    GetRMS(rms, rmsTotal);
    BlockIndices maxlocs = MaxLocs();
    o << "RMS: ";
    for (int i = 0; i < mDOF; ++i)
    {
        o << rms[i] << "(" << maxlocs[i] << ") ";
    }
    o << ", total = " << rmsTotal << "(" << maxlocs[mDOF] << ")";
    delete[] rms;
    return o;
}

Residual
Residual::operator + (const Residual& r) const
{
    Residual result(*this);

    for (int i = 0; i < mDOF; ++i)
    {
        result.mSquare[i] = result.mSquare[i] + r.mSquare[i];
    }
    for (int i = 0; i < mDOF + 1; ++i)
    {
        if (result.mMaxSquare[i] < r.mMaxSquare[i])
        {
            result.mMaxSquare[i] = r.mMaxSquare[i];
            result.mMaxLocs[i] = r.mMaxLocs[i];
        }
    }

    result.mNumPoints += r.mNumPoints;

    return result;
}

