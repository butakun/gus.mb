// $Id: MPIFuncs.cpp 191 2012-01-20 02:38:43Z kato $

#include "MPIFuncs.h"
#include <mpi.h>
#include <cassert>

void Init()
{
    int argc = 1;
    char** argv;

    int err;
    err = MPI_Init(&argc, &argv);
    assert(err == MPI_SUCCESS);
}

void Finalize()
{
    int err;
    err = MPI_Finalize();
    assert(err == MPI_SUCCESS);
}

int Rank()
{
    int err;
    int rank;
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(err == MPI_SUCCESS);
    return rank;
}

