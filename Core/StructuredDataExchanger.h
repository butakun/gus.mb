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
// $Id: StructuredDataExchanger.h 270 2013-02-08 09:26:36Z kato $
#ifndef INCLUDED_STRUCTURED_DATA_EXCHANGER_H__
#define INCLUDED_STRUCTURED_DATA_EXCHANGER_H__

#include "DataExchanger.h"
#include "Block.h"
#include "PatchExchanger.h"
#include <list>
#include <map>

class StructuredDataExchanger : public DataExchanger
{
public:
    typedef std::list<PatchExchanger*> PatchExchangers;
    typedef std::map<int, Structured<double>*> BlockDataMap;

    StructuredDataExchanger(const Block& block, Structured<double>& data);
    virtual ~StructuredDataExchanger();

    virtual void Start();
    virtual void Finish();

    Structured<double>& Data() const { return mData; }

    // used by PatchExchanger
    Structured<double>* BlockData(int block);

protected:

private:
    static BlockDataMap mBlockDataMap;

    const Block& mBlock;
    Structured<double>& mData;
    PatchExchangers mPatchExchangers;
};

#endif // INCLUDED_STRUCTURED_DATA_EXCHANGER_H__

