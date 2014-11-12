// $Id: CGNSWriter.h 294 2013-08-23 14:11:39Z kato $
#ifndef INCLUDE_CGNS_WRITER_H__
#define INCLUDE_CGNS_WRITER_H__

#include "Structured.h"
#include "Block.h"
#include "Physics.h"
#include <string>

class CGNSStructure;

class CGNSWriter
{
public:
    CGNSWriter(const char* filename, bool modify = true);
    ~CGNSWriter();

    void WriteStructure(const CGNSStructure& s);
    void WriteFlowSolution(int zone, const Block& block, const Structured<double>& U, const Physics& phys);
    void WriteTurbulenceSolution(int zone, const Block& block, const Structured<double>& U, const Structured<double>& UT, const Physics& phys, const std::string& model);

    int WriteBase();
    int WriteZone(int B, const IndexRange& meshRange, const char* name = NULL);

protected:

private:
    int fn;
};

#endif // INCLUDE_CGNS_WRITER_H__

