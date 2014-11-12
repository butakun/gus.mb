// $Id: PLOT3DMeshReader.h 4 2010-03-06 14:10:00Z kato $
#ifndef INCLDUED_HIRO_PLOT3D_MESH_READER_H__
#define INCLDUED_HIRO_PLOT3D_MESH_READER_H__

#include "MeshReader.h"
#include <string>

class PLOT3DMeshReader : public MeshReader
{
public:
    PLOT3DMeshReader(const char* filename);
    virtual ~PLOT3DMeshReader();

    virtual void Read(Structured<double>& xyz, IndexRange& range);

protected:

private:
    std::string mFileName;
};

#endif // INCLDUED_HIRO_PLOT3D_MESH_READER_H__

