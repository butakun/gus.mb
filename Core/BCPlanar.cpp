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
// $Id: BCPlanar.cpp 231 2012-05-11 01:11:28Z kato $

#include "BCPlanar.h"
#include "IndexUtils.h"
#include <cassert>

IndexRange
BCPlanar::RindRange() const
{
    return IndexUtils::FromMeshRangeToCellRange(MeshRange(), PatchDirection(), 1); // FIXME: only single ghost layer?
}

void
BCPlanar::SetMaskWithAValue(Structured<int>& mask, int value) const
{
    IndexRange rindRange = RindRange();
    for (IndexIterator itor(rindRange); !itor.IsEnd(); itor.Advance())
    {
        IndexIJK ijk = itor.Index();
        *mask(ijk) = value;
    }
}

