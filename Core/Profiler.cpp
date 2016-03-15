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
// $Id: Profiler.cpp 189 2012-01-19 07:40:08Z kato $

#include <mpi.h>
#include "Profiler.h"
#include <cstdlib>
#include <cassert>

Profiler* Profiler::mInstance = NULL;

Profiler*
Profiler::GetInstance()
{
    if (mInstance == NULL)
    {
        mInstance = new Profiler();
    }

    return mInstance;
}

Profiler::Profiler()
{
    mTimeStamp = MPI_Wtime();
}

Profiler::~Profiler()
{
}

double
Profiler::CheckPoint()
{
    double t0 = mTimeStamp;
    mTimeStamp = MPI_Wtime();
    return mTimeStamp - t0;
}

