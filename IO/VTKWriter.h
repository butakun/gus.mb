/*
    gus.mb, an open source flow solver.
    Copyright (C) 2016 Hiromasa Kato <hiromasa at gmail.com>

    This file is part of gus.mb.

    gus.mb is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    gus.mb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
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

