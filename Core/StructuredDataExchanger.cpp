// $Id: StructuredDataExchanger.cpp 180 2012-01-13 08:48:04Z kato $

#include <mpi.h> // FIXME
#include "StructuredDataExchanger.h"
#include "Communicator.h"
#include "Connectivity1to1.h"
#include "PatchExchanger1to1.h"
#include <vector>

StructuredDataExchanger::BlockDataMap StructuredDataExchanger::mBlockDataMap;

StructuredDataExchanger::StructuredDataExchanger(const Block& block, Structured<double>& data)
:   mBlock(block), mData(data)
{
    const BCs& bcs = mBlock.GetBCs();
    for (BCs::const_iterator i = bcs.begin(); i != bcs.end(); ++i)
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

