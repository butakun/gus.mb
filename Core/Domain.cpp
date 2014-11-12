// $Id: Domain.cpp 87 2011-01-24 12:43:58Z kato $

#include "Domain.h"
#include "Block.h"
#include "Roster.h"
#include "Communicator.h"
#include "BCViscousWall.h"
#include <algorithm>
#include <vector>

void
Domain::RegisterLocalBlock(Block* block)
{
    mLocalBlocks.push_back(block);
    int rank = Communicator::GetInstance()->MyRank();
    Roster::GetInstance()->RegisterBlock(rank, block->ID());
}

void
Domain::RegisterRemoteBlock(int rank, int blockID)
{
    Roster::GetInstance()->RegisterBlock(rank, blockID);
}

Block*
Domain::FindLocalBlock(int zone)
{
    Blocks::iterator i = std::find_if(mLocalBlocks.begin(), mLocalBlocks.end(),
        std::bind2nd(std::mem_fun(&Block::IsZone), zone));
    assert(i != mLocalBlocks.end());
    return *i;
}

#include <mpi.h> // FIXME: not mpi dependency here

typedef struct { const Block* block; IndexRange meshRange; } wall_s;

void
Domain::GatherWalls()
{
    std::vector<wall_s> bcwalls;

    for (Blocks::const_iterator i = mLocalBlocks.begin(); i != mLocalBlocks.end(); ++i)
    {
        const Block* block = *i;
        const BCs& bcs = block->GetBCs();
        for (BCs::const_iterator j = bcs.begin(); j != bcs.end(); ++j)
        {
            const BC* bc = *j;
            const BCViscousWall* bcwall = dynamic_cast<const BCViscousWall*>(bc);
            if (bcwall == NULL)
            {
                continue;
            }
            std::cout << "Domain::GatherWalls: wall = " << bcwall->MeshRange() << std::endl;
            wall_s ws = { block, bcwall->MeshRange() };
            bcwalls.push_back(ws);
        }
    }

    int err;
    int sendbuf = bcwalls.size();
    int* recvbuf = new int[Communicator::GetInstance()->Size()];
    err = MPI_Allgather(&sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);

    delete[] recvbuf;
}

