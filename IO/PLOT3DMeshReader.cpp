// $Id: PLOT3DMeshReader.cpp 4 2010-03-06 14:10:00Z kato $

#include "PLOT3DMeshReader.h"
#include <cassert>

PLOT3DMeshReader::PLOT3DMeshReader(const char* filename)
:   mFileName(filename)
{
}

PLOT3DMeshReader::~PLOT3DMeshReader()
{
}

void
PLOT3DMeshReader::Read(Structured<double>& xyz, IndexRange& range)
{
    assert(false);
}

