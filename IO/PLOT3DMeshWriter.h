// $Id: PLOT3DMeshWriter.h 76 2010-10-06 16:50:49Z kato $
#ifndef INCLUDED_PLOT3D_MESH_WRITER_H__
#define INCLUDED_PLOT3D_MESH_WRITER_H__

#include <vector>

class Block;

class PLOT3DMeshWriter
{
public:
    PLOT3DMeshWriter() {}

    void AddBlock(const Block& block);
    void Write(const char* filename) const;

    void Write(const char* filename, const Block& block) const;

protected:

private:
    std::vector<const Block*> mBlocks;
};

#endif // INCLUDED_PLOT3D_MESH_WRITER_H__

