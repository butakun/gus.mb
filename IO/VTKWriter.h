// $Id: VTKWriter.h 250 2012-07-20 01:52:56Z kato $
#ifndef INCLUDED_VTK_WRITER_H__
#define INCLUDED_VTK_WRITER_H__

#include "Structured.h"
#include <string>
#include <vector>

class VTKWriter
{
public:
    enum Format { LEGACY, XML };
    enum Type { SCALAR, VECTOR };

    VTKWriter(const char* filename, Format format = LEGACY);
    ~VTKWriter();

    void AddMesh(const Structured<double>& xyz);
    void AddData(const Structured<double>& data, int offset, const char* name, Type type, double scale = 1.0);

    void Write();
    void Finalize(const char* filename, const std::vector<std::string>& filenames);

protected:
    void WriteLegacy();
    void WriteScalar(std::ostream& o, const char* name, const Structured<double>& data, int offset, double scale);
    void WriteVector(std::ostream& o, const char* name, const Structured<double>& data, int offset, double scale);

    void WriteXML();
    void WriteScalarXML(std::ostream& o, const char* name, const Structured<double>& data, int offset, double scale);
    void WriteVectorXML(std::ostream& o, const char* name, const Structured<double>& data, int offset, double scale);

private:
    typedef struct {
        std::string name;
        int offset;
        double scale;
        Type type;
        const Structured<double>* data;
    } DataInfo;
    typedef std::vector<DataInfo> Data;

    std::string mFileName;
    Format mFormat;
    const Structured<double>* mXYZ;
    Data mData;
};

#endif // INCLUDED_VTK_WRITER_H__

