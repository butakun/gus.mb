// $Id: MeshReader.h 4 2010-03-06 14:10:00Z kato $
#ifndef INCLUDED_HIRO_MESH_READER_H__
#define INCLUDED_HIRO_MESH_READER_H__

#include "Structured.h"

class MeshReader
{
public:
    virtual ~MeshReader() {};

    virtual void Read(Structured<double>& xyz, IndexRange& range) = 0;

protected:

private:
};

#endif // INCLUDED_HIRO_MESH_READER_H__

