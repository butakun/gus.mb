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

