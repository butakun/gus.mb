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

#include "BCPlanar.h"
#include "IndexUtils.h"
#include <cassert>

BCPlanar::BCPlanar(const IndexRange& meshRange, Direction direction)
:   BC(meshRange), mDirection(direction),
    mDGhost(0, 0, 0), mDInterior(0, 0, 0), mDeltaInterior(0, 0, 0),
    mSnDirection(1.0)
{
    switch (PatchDirection())
    {
    case I:
        mDInterior.I = 1;
        mDeltaInterior.I = 1;
        mSnDirection = 1.0;
        break;
    case INEG:
        mDGhost.I = 1;
        mDeltaInterior.I = -1;
        mSnDirection = -1.0;
        break;
    case J:
        mDInterior.J = 1;
        mDeltaInterior.J = 1;
        mSnDirection = 1.0;
        break;
    case JNEG:
        mDGhost.J = 1;
        mDeltaInterior.J = -1;
        mSnDirection = -1.0;
        break;
    case K:
        mDInterior.K = 1;
        mDeltaInterior.K = 1;
        mSnDirection = 1.0;
        break;
    case KNEG:
        mDGhost.K = 1;
        mDeltaInterior.K = -1;
        mSnDirection = -1.0;
        break;
    default:
        assert(false);
    }
}

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

