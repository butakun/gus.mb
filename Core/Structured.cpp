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
// $Id: Structured.cpp 14 2010-04-14 10:56:36Z kato $

#include "Structured.h"

IndexRange
IndexRange::Canonical() const
{
    int imin, jmin, kmin, imax, jmax, kmax;
    imin = std::min(Start.I, End.I);
    imax = std::max(Start.I, End.I);
    jmin = std::min(Start.J, End.J);
    jmax = std::max(Start.J, End.J);
    kmin = std::min(Start.K, End.K);
    kmax = std::max(Start.K, End.K);

    return IndexRange(imin, jmin, kmin, imax, jmax, kmax);
}

