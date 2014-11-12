// $Id: CGNSMeshReader.h 76 2010-10-06 16:50:49Z kato $
#ifndef INCLUDED_CGNS_MESH_READER_H__
#define INCLUDED_CGNS_MESH_READER_H__

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "DomainInfo.h"
#include "Structured.h"

class Block;
class BC;

class CGNSMeshReader
{
public:
    CGNSMeshReader(const char* filename);
    virtual ~CGNSMeshReader();

    const DomainInfo& Domain() const { return mDomain; }
    void ReadMesh(Structured<double>& XYZ, int Z) const;
    void ReadZone(Block** block, int Z);
    void ReadBC(BC** bc, int Z, int BC);
    void ReadPartialMesh(Structured<double>& XYZ, int Z, const IndexRange& range) const;

protected:
    void ReadBlockInfo();
    void GetZone(int& fn, int B, int Z) const;
    void Close(int fn) const;

private:
    std::string mFileName;
    DomainInfo mDomain;
};

#endif // INCLUDED_CGNS_MESH_READER_H__

