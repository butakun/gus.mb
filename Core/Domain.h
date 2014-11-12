// $Id: Domain.h 87 2011-01-24 12:43:58Z kato $
#ifndef INCLUDED_DOMAIN_H__
#define INCLUDED_DOMAIN_H__

#include "DomainInfo.h"
#include <vector>

class Block;
class Wall;

typedef std::vector<Block*> Blocks;
typedef std::vector<Wall*> Walls;

class Domain
{
public:
    Domain(const DomainInfo& info) : mInfo(info) {}
    virtual ~Domain() {}

    const DomainInfo& GetDomainInfo() const { return mInfo; }
    void RegisterLocalBlock(Block* block); // { mLocalBlocks.push_back(block); }
    void RegisterRemoteBlock(int rank, int blockID);
    const Blocks& GetLocalBlocks() const { return mLocalBlocks; }
    Block* FindLocalBlock(int zone);

    void GatherWalls();
    const Walls& GetWalls() const { return mWalls; }

protected:

private:
    DomainInfo mInfo;
    Blocks mLocalBlocks;
    Walls mWalls;
};

#endif // INCLUDED_DOMAIN_H__

