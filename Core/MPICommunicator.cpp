// $Id: MPICommunicator.cpp 225 2012-05-01 23:18:07Z kato $

#include <mpi.h>
#include "MPICommunicator.h"
#include <sstream>
#include <cassert>

MPICommunicator::MPICommunicator(int* argc, char*** argv)
:   mNextTag(1234)
{
    int err;
    err = MPI_Init(argc, argv);
    assert(err == MPI_SUCCESS);

    std::ostringstream oss;
    oss << "log." << MyRank();
    mConsole.open(oss.str().c_str());

    mConsole << "Communicator size = " << Size() << std::endl;
}

MPICommunicator::~MPICommunicator()
{
    mConsole.close();
}

int
MPICommunicator::Size() const
{
    int err;
    int size;
    err = MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(err == MPI_SUCCESS);
    return size;
}

int
MPICommunicator::MyRank() const
{
    int err;
    int rank;
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(err == MPI_SUCCESS);
    return rank;
}

//MPICommunicator::

void
MPICommunicator::Barrier() const
{
    int err;
    err = MPI_Barrier(MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
}

void
MPICommunicator::Finalize() const
{
    int err;
    err = MPI_Finalize();
    assert(err == MPI_SUCCESS);
}

std::ostream&
MPICommunicator::Console()
{
    return mConsole;
}

int
MPICommunicator::ReserveTag(void* key)
{
    TagMap::const_iterator i = mTagMap.find(key);
    if (i == mTagMap.end())
    {
        // A new tag
        int tag = mNextTag++;
        mTagMap[key] = tag;
        return tag;
    }
    return i->second;
}

void
MPICommunicator::Broadcast(int* buffer, int count, int root) const
{
    int err;
    err = MPI_Bcast(buffer, count, MPI_INT, root, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
}

bool
MPICommunicator::Any(bool flag) const
{
    int err;
    int send = flag, recv;
    err = MPI_Allreduce(&send, &recv, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);

    return recv != 0;
}

int
MPICommunicator::ReduceSum(int value, int root) const
{
    int err;
    int recv;
    err = MPI_Reduce(&value, &recv, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);

    return recv;
}

int
MPICommunicator::AllReduceSum(int value) const
{
    int err;
    int recv;
    err = MPI_Allreduce(&value, &recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);

    return recv;
}

