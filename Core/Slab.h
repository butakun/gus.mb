// $Id: Slab.h 41 2010-07-06 10:48:07Z kato $
#ifndef INCLUDED_SLAB_H__
#define INCLUDED_SLAB_H__

#include "Structured.h"
#include <string>
#include <map>

class Slab
{
public:
    Slab(int blockID, const char* name);
    virtual ~Slab();

    int BlockID() const { return mBlockID; }
    const std::string& Name() const { return mName; }

    Structured<double>& AddData(const char* name, int dof, const IndexRange& range);
    Structured<double>& GetData(const char* name);

protected:

private:
    typedef std::map<std::string, Structured<double> > DataMap;

    int mBlockID;
    std::string mName;
    DataMap mDataMap;
};

#endif // INCLUDED_SLAB_H__

