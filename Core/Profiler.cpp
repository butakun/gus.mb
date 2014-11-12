// $Id: Profiler.cpp 189 2012-01-19 07:40:08Z kato $

#include <mpi.h>
#include "Profiler.h"
#include <cstdlib>
#include <cassert>

Profiler* Profiler::mInstance = NULL;

Profiler*
Profiler::GetInstance()
{
    if (mInstance == NULL)
    {
        mInstance = new Profiler();
    }

    return mInstance;
}

Profiler::Profiler()
{
    mTimeStamp = MPI_Wtime();
}

Profiler::~Profiler()
{
}

double
Profiler::CheckPoint()
{
    double t0 = mTimeStamp;
    mTimeStamp = MPI_Wtime();
    return mTimeStamp - t0;
}

