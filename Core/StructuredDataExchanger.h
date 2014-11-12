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

