// $Id: PLOT3DSolutionWriter.h 78 2010-12-09 13:44:56Z kato $
#ifndef INCLUDED_PLOT3D_SOLUTION_WRITER_H__
#define INCLUDED_PLOT3D_SOLUTION_WRITER_H__

#include "Structured.h"
#include <string>

class PLOT3DSolutionWriter
{
public:
    PLOT3DSolutionWriter(const char* name);

    void Write(const Structured<double>& U, const IndexRange& meshRange);

protected:

private:
    std::string mFileName;
};

#endif // INCLUDED_PLOT3D_SOLUTION_WRITER_H__

