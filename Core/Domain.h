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

