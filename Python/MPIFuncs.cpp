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
// $Id: MPIFuncs.cpp 191 2012-01-20 02:38:43Z kato $

#include "MPIFuncs.h"
#include <mpi.h>
#include <cassert>

void Init()
{
    int argc = 1;
    char** argv;

    int err;
    err = MPI_Init(&argc, &argv);
    assert(err == MPI_SUCCESS);
}

void Finalize()
{
    int err;
    err = MPI_Finalize();
    assert(err == MPI_SUCCESS);
}

int Rank()
{
    int err;
    int rank;
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(err == MPI_SUCCESS);
    return rank;
}

