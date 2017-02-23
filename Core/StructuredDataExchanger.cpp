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

#include <mpi.h> // FIXME
#include "StructuredDataExchanger.h"
#include "Communicator.h"
#include "Connectivity1to1.h"
#include "PatchExchanger1to1.h"
#include "Roster.h"
#include <vector>

StructuredDataExchanger::BlockDataMap StructuredDataExchanger::mBlockDataMap;

StructuredDataExchanger::StructuredDataExchanger(const Block& block, const Model* model, Structured<double>& data)
:   mBlock(block), mData(data)
{
    const std::vector<BC*>& bcs = Roster::GetInstance()->GetBCs(block.ID(), model);
    for (std::vector<BC*>::const_iterator i = bcs.begin(); i != bcs.end(); ++i)
    {
        Connectivity* conn = dynamic_cast<Connectivity*>(*i);
        if (conn == NULL)
            continue;
        Connectivity1to1* conn1to1 = dynamic_cast<Connectivity1to1*>(conn);
        if (conn1to1 == NULL)
            assert(false);
        PatchExchanger1to1* patchEx = new PatchExchanger1to1(this, *conn1to1, mData);
        mPatchExchangers.push_back(patchEx);
    }
}

StructuredDataExchanger::~StructuredDataExchanger()
{
    for (PatchExchangers::iterator i = mPatchExchangers.begin();
         i != mPatchExchangers.end(); ++i)
    {
        delete *i;
    }
}

void
StructuredDataExchanger::Start()
{
    //std::cout << "Registering exchange data " << &mData << " of Block " << mBlock.ID() << std::endl;
    mBlockDataMap[mBlock.ID()] = &mData;

    for (PatchExchangers::iterator i = mPatchExchangers.begin();
         i != mPatchExchangers.end(); ++i)
    {
        (*i)->Start();
    }
}

void
StructuredDataExchanger::Finish()
{
    for (PatchExchangers::iterator i = mPatchExchangers.begin();
         i != mPatchExchangers.end(); ++i)
    {
        (*i)->Finish();
    }
}

Structured<double>*
StructuredDataExchanger::BlockData(int block)
{
    BlockDataMap::iterator i = mBlockDataMap.find(block);

    Structured<double>* data = NULL;
    if (i != mBlockDataMap.end())
    {
        data = i->second;
    }
    return data;
}

