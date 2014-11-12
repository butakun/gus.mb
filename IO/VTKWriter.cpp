// $Id: VTKWriter.cpp 250 2012-07-20 01:52:56Z kato $

#include "VTKWriter.h"
#include "IndexUtils.h"
#include <fstream>

VTKWriter::VTKWriter(const char* filename, Format format)
:   mFileName(filename), mFormat(format), mXYZ(NULL)
{
}

VTKWriter::~VTKWriter()
{
}


void
VTKWriter::AddMesh(const Structured<double>& xyz)
{
    assert(mXYZ == NULL);
    mXYZ = &xyz;
}

void
VTKWriter::AddData(const Structured<double>& data, int offset, const char* name, Type type, double scale)
{
    DataInfo info;
    info.name = name;
    info.type = type;
    info.offset = offset;
    info.scale = scale;
    info.data = &data;
    mData.push_back(info);
}

void
VTKWriter::Write()
{
    if (mFormat == LEGACY)
        WriteLegacy();
    else if (mFormat == XML)
        WriteXML();
    else
        assert(false);
}

void
VTKWriter::WriteLegacy()
{
    assert(mXYZ != NULL);

    const Structured<double>& xyz = *mXYZ;
    IndexRange mr = xyz.GetRange();
    IndexIJK ms = mr.Shape();

    std::ofstream o(mFileName.c_str());

    o << "# vtk DataFile Version 2.0" << std::endl;
    o << "Noname" << std::endl;
    o << "ASCII" << std::endl;
    o << "DATASET STRUCTURED_GRID" << std::endl;
    o << "DIMENSIONS " << ms.I << ' ' << ms.J << ' ' << ms.K << std::endl;
    o << "POINTS " << mr.Count() << " float" << std::endl;
    for (int k = mr.Start.K; k <= mr.End.K; ++k)
    {
        for (int j = mr.Start.J; j <= mr.End.J; ++j)
        {
            for (int i = mr.Start.I; i <= mr.End.I; ++i)
            {
                double* p = xyz(i, j, k);
                o << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
            }
        }
    }

    if (!mData.empty())
    {
        IndexRange mr = mXYZ->GetRange();
        IndexRange cr(mr.Start.I + 1, mr.Start.J + 1, mr.Start.K + 1, mr.End.I, mr.End.J, mr.End.K);
        o << "CELL_DATA " << cr.Count() << std::endl;
    }
    for (Data::const_iterator i = mData.begin(); i != mData.end(); ++i)
    {
        const DataInfo& di = *i;
        if (di.type == SCALAR)
        {
            WriteScalar(o, di.name.c_str(), *di.data, di.offset, di.scale);
        }
        else if (di.type == VECTOR)
        {
            WriteVector(o, di.name.c_str(), *di.data, di.offset, di.scale);
        }
        else
        {
            assert(false);
        }
    }
}

void
VTKWriter::WriteScalar(std::ostream& o, const char* name, const Structured<double>& data, int offset, double scale)
{
    IndexRange mr = mXYZ->GetRange();
    IndexRange cr(mr.Start.I + 1, mr.Start.J + 1, mr.Start.K + 1, mr.End.I, mr.End.J, mr.End.K);
    o << "SCALARS " << name << " float 1" << std::endl;
    o << "LOOKUP_TABLE default" << std::endl;
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* d = data(i, j, k);
                o << d[offset] * scale << std::endl;
            }
        }
    }
}

void
VTKWriter::WriteVector(std::ostream& o, const char* name, const Structured<double>& data, int offset, double scale)
{
    IndexRange mr = mXYZ->GetRange();
    IndexRange cr(mr.Start.I + 1, mr.Start.J + 1, mr.Start.K + 1, mr.End.I, mr.End.J, mr.End.K);
    o << "VECTORS " << name << " float" << std::endl;
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* d = &data(i, j, k)[offset];
                o << d[0] * scale << " " << d[1] * scale << " " << d[2] * scale << std::endl;
            }
        }
    }
}

void
VTKWriter::WriteXML()
{
    assert(mXYZ != NULL);

    const Structured<double>& xyz = *mXYZ;
    IndexRange mr = xyz.GetRange();
    IndexIJK ms = mr.Shape();

    std::ofstream o(mFileName.c_str());

    o << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    o << "  <StructuredGrid WholeExtent=\"" << mr.Start.I << ' ' << mr.End.I << ' ' << mr.Start.J << ' ' << mr.End.J << ' ' << mr.Start.K << ' ' << mr.End.K << "\">" << std::endl;
    o << "    <Piece Extent=\"" << mr.Start.I << ' ' << mr.End.I << ' ' << mr.Start.J << ' ' << mr.End.J << ' ' << mr.Start.K << ' ' << mr.End.K << "\">" << std::endl;

    // Grid points
    o << "      <Points>" << std::endl;
    o << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    int two = 0;
    for (int k = mr.Start.K; k <= mr.End.K; ++k)
    {
        for (int j = mr.Start.J; j <= mr.End.J; ++j)
        {
            for (int i = mr.Start.I; i <= mr.End.I; ++i)
            {
                double* p = xyz(i, j, k);
                o << p[0] << ' ' << p[1] << ' ' << p[2];
                ++two;
                if (two % 2 == 0)
                    o << std::endl;
                else
                    o << ' ';
            }
        }
    }
    if (two % 2 != 0)
        o << std::endl;
    o << "        </DataArray>" << std::endl;
    o << "      </Points>" << std::endl;

    if (!mData.empty())
    {
        IndexRange mr = mXYZ->GetRange();
        IndexRange cr(mr.Start.I + 1, mr.Start.J + 1, mr.Start.K + 1, mr.End.I, mr.End.J, mr.End.K);
        o << "      <CellData>" << std::endl;

        for (Data::const_iterator i = mData.begin(); i != mData.end(); ++i)
        {
            const DataInfo& di = *i;
            if (di.type == SCALAR)
            {
                WriteScalarXML(o, di.name.c_str(), *di.data, di.offset, di.scale);
            }
            else if (di.type == VECTOR)
            {
                WriteVectorXML(o, di.name.c_str(), *di.data, di.offset, di.scale);
            }
            else
            {
                assert(false);
            }
        }

        o << "      </CellData>" << std::endl;
    }
    o << "    </Piece>>" << std::endl;
    o << "  </StructuredGrid>" << std::endl;
    o << "</VTKFile>" << std::endl;
}

void
VTKWriter::WriteScalarXML(std::ostream& o, const char* name, const Structured<double>& data, int offset, double scale)
{
    IndexRange mr = mXYZ->GetRange();
    IndexRange cr(mr.Start.I + 1, mr.Start.J + 1, mr.Start.K + 1, mr.End.I, mr.End.J, mr.End.K);
    o << "<DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">" << std::endl;
    int six = 0;
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* d = data(i, j, k);
                o << d[offset] * scale;
                six++;
                if (six % 6 == 0)
                    o << std::endl;
                else
                    o << ' ';
            }
        }
    }
    o << "</DataArray>" << std::endl;
}

void
VTKWriter::WriteVectorXML(std::ostream& o, const char* name, const Structured<double>& data, int offset, double scale)
{
    IndexRange mr = mXYZ->GetRange();
    IndexRange cr(mr.Start.I + 1, mr.Start.J + 1, mr.Start.K + 1, mr.End.I, mr.End.J, mr.End.K);
    o << "<DataArray type=\"Float32\" Name=\"" << name << "\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    int two = 0;
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* d = &data(i, j, k)[offset];
                o << d[0] * scale << " " << d[1] * scale << " " << d[2] * scale ;
                two++;
                if (two % 2 == 0)
                    o << std::endl;
                else
                    o << ' ';
            }
        }
    }
    o << "</DataArray>" << std::endl;
}

void
VTKWriter::Finalize(const char* filename, const std::vector<std::string>& filenames)
{
    assert(mFormat == XML);

    std::ofstream o(filename);

    o << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl;
    o << "<vtkMultiBlockDataSet>" << std::endl;
    for (size_t i = 0; i < filenames.size(); ++i)
    {
        o << "<DataSet index=\"" << i << "\" file=\"" << filenames[i] << "\"/>" << std::endl;
    }
    o << "</vtkMultiBlockDataSet>" << std::endl;
    o << "</VTKFile>" << std::endl;
}

