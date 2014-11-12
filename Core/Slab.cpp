// $Id: Slab.cpp 41 2010-07-06 10:48:07Z kato $

#include "Slab.h"
#include <cassert>

Slab::Slab(int blockID, const char* name)
:   mBlockID(blockID), mName(name)
{
}

Slab::~Slab()
{
    for (DataMap::iterator i = mDataMap.begin(); i != mDataMap.end(); ++i)
    {
        Structured<double>& d = i->second;
        delete[] d.Data;
    }
}

Structured<double>&
Slab::AddData(const char* name, int dof, const IndexRange& range)
{
    DataMap::const_iterator i = mDataMap.find(name);
    assert(i == mDataMap.end());

    Structured<double> d(dof, range);
    mDataMap[name] = d;

    return mDataMap[name];
}

Structured<double>&
Slab::GetData(const char* name)
{
    DataMap::iterator i = mDataMap.find(name);
    assert(i != mDataMap.end());

    return i->second;
}

