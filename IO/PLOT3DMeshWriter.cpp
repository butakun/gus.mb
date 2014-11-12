// $Id: PLOT3DMeshWriter.cpp 89 2011-01-28 14:16:26Z kato $

#include "PLOT3DMeshWriter.h"
#include "Block.h"
#include <fstream>

void
PLOT3DMeshWriter::AddBlock(const Block& block)
{
    mBlocks.push_back(&block);
}

void WriteBlock(std::ostream& o, const Block& block)
{
    IndexRange mr = block.MeshRange();
    const Structured<double>& xyz = block.XYZ();

    for (int l = 0; l < 3; ++l)
    {
        for (int k = mr.Start.K; k <= mr.End.K; ++k)
        {
            for (int j = mr.Start.J; j <= mr.End.J; ++j)
            {
                for (int i = mr.Start.I; i <= mr.End.I; ++i)
                {
                    o << xyz(i, j, k)[l] << std::endl;
                }
            }
        }
    }
}

void
PLOT3DMeshWriter::Write(const char* filename) const
{
    std::ofstream f(filename);

    f << mBlocks.size() << std::endl;
    for (size_t i = 0; i < mBlocks.size(); ++i)
    {
        const Block& block = *mBlocks[i];
        IndexIJK shape = block.MeshRange().Shape();
        f << shape.I << '\t' << shape.J << '\t' << shape.K << std::endl;
    }
    for (size_t i = 0; i < mBlocks.size(); ++i)
    {
        const Block& block = *mBlocks[i];
        WriteBlock(f, block);
    }

    f.close();
}

void
PLOT3DMeshWriter::Write(const char* filename, const Block& block) const
{
    std::ofstream f(filename);

    IndexIJK shape = block.MeshRange().Shape();
    f << shape.I << '\t' << shape.J << '\t' << shape.K << std::endl;
    WriteBlock(f, block);

    f.close();
}

