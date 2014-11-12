// $Id$
#ifndef INCLUDED_PLOT3D_FUNCTION_WRITER_H__
#define INCLUDED_PLOT3D_FUNCTION_WRITER_H__

#include <string>
#include <iostream>

class PLOT3DFunctionWriter
{
public:
    PLOT3DFunctionWriter(const char* name);

    void Write(const Structured<double>& f, double* (*func)(double *));

protected:

private:
    std::string mFileName;
    std::ostream& mOut;
};

#endif // INCLUDED_PLOT3D_FUNCTION_WRITER_H__

