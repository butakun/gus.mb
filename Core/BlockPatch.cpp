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

#include "Roster.h"
#include "BlockPatch.h"
#include "IndexUtils.h"
#include <cassert>

BlockPatch
BlockPatch::New(int blockID, const IndexRange& meshRange, int uniqueID)
{
    BlockPatch bp(blockID, meshRange, uniqueID);
    BlockPatch bp2 = Roster::GetInstance()->RegisterBlockPatch(bp);
    return bp2;
}

void
BlockPatch::InitOrientation()
{
    assert(mMeshRange.IsCanonical());
    IndexRange bmr = Roster::GetInstance()->GetBlock(mBlockID)->MeshRange();
    mDir = IndexUtils::PatchDirection(bmr, mMeshRange);

    switch (mDir)
    {
    case I:
        mI1 = IndexIJK(0, 1, 0);
        mI2 = IndexIJK(0, 0, 1);
        mI3Ghost    = IndexIJK(0, 0, 0);
        mI3Interior = IndexIJK(1, 0, 0);
        mDI3Ghost = IndexIJK(-1, 0, 0);
        break;
    case INEG:
        mI1 = IndexIJK(0, 1, 0);
        mI2 = IndexIJK(0, 0, 1);
        mI3Ghost    = IndexIJK(1, 0, 0);
        mI3Interior = IndexIJK(0, 0, 0);
        mDI3Ghost = IndexIJK(1, 0, 0);
        break;
    case J:
        mI1 = IndexIJK(1, 0, 0);
        mI2 = IndexIJK(0, 0, 1);
        mI3Ghost    = IndexIJK(0, 0, 0);
        mI3Interior = IndexIJK(0, 1, 0);
        mDI3Ghost = IndexIJK(0, -1, 0);
        break;
    case JNEG:
        mI1 = IndexIJK(1, 0, 0);
        mI2 = IndexIJK(0, 0, 1);
        mI3Ghost    = IndexIJK(0, 1, 0);
        mI3Interior = IndexIJK(0, 0, 0);
        mDI3Ghost = IndexIJK(0, 1, 0);
        break;
    case K:
        mI1 = IndexIJK(1, 0, 0);
        mI2 = IndexIJK(0, 1, 0);
        mI3Ghost    = IndexIJK(0, 0, 0);
        mI3Interior = IndexIJK(0, 0, 1);
        mDI3Ghost = IndexIJK(0, 0, -1);
        break;
    case KNEG:
        mI1 = IndexIJK(1, 0, 0);
        mI2 = IndexIJK(0, 1, 0);
        mI3Ghost    = IndexIJK(0, 0, 1);
        mI3Interior = IndexIJK(0, 0, 0);
        mDI3Ghost = IndexIJK(0, 0, 1);
        break;
    default:
        assert(false);
    }
}

IndexRange
BlockPatch::CellFaceRange() const
{
    IndexRange cfr = mMeshRange;
    switch (mDir)
    {
    case I:
    case INEG:
        cfr.Start += IndexIJK(0, 1, 1);
        break;
    case J:
    case JNEG:
        cfr.Start += IndexIJK(1, 0, 1);
        break;
    case K:
    case KNEG:
        cfr.Start += IndexIJK(1, 1, 0);
        break;
    default:
        assert(false);
    }

    return cfr;
}

void
BlockPatch::GhostCellRange(IndexRange& gcr, IndexIJK& i1, IndexIJK& i2, IndexIJK& di3) const
{
    assert(mMeshRange.IsCanonical());
    IndexRange bmr = Roster::GetInstance()->GetBlock(mBlockID)->MeshRange();
    Direction dir = IndexUtils::PatchDirection(bmr, mMeshRange);

    gcr = mMeshRange;
    switch (dir)
    {
    case I:
        gcr.Start += IndexIJK(0, 1, 1);
        i1 = IndexIJK(0, 1, 0);
        i2 = IndexIJK(0, 0, 1);
        di3 = IndexIJK(-1, 0, 0);
        break;
    case INEG:
        gcr.Start += IndexIJK(1, 1, 1);
        gcr.End += IndexIJK(1, 0, 0);
        i1 = IndexIJK(0, 1, 0);
        i2 = IndexIJK(0, 0, 1);
        di3 = IndexIJK(1, 0, 0);
        break;
    case J:
        gcr.Start += IndexIJK(1, 0, 1);
        i1 = IndexIJK(1, 0, 0);
        i2 = IndexIJK(0, 0, 1);
        di3 = IndexIJK(0, -1, 0);
        break;
    case JNEG:
        gcr.Start += IndexIJK(1, 1, 1);
        gcr.End += IndexIJK(0, 1, 0);
        i1 = IndexIJK(1, 0, 0);
        i2 = IndexIJK(0, 0, 1);
        di3 = IndexIJK(0, 1, 0);
        break;
    case K:
        gcr.Start += IndexIJK(1, 1, 0);
        i1 = IndexIJK(1, 0, 0);
        i2 = IndexIJK(0, 1, 0);
        di3 = IndexIJK(0, 0, -1);
        break;
    case KNEG:
        gcr.Start += IndexIJK(1, 1, 1);
        gcr.End += IndexIJK(0, 0, 1);
        i1 = IndexIJK(1, 0, 0);
        i2 = IndexIJK(0, 1, 0);
        di3 = IndexIJK(0, 0, 1);
        break;
    default:
        assert(false);
    }
}

IndexRange
BlockPatch::MeshRangeToCellRange(const IndexRange& mr) const
{
    assert(mr.IsCanonical());

    IndexRange cr = mr;
    IndexIJK ms = mr.Shape();

    IndexIJK d(0, 0, 0);
    IndexIJK ds(0, 0, 0);

    // FIXME: should use IndexUtils::PatchDirection (so it won't depend on Start.I being zero, etc.)
    if (ms.I == 1)
    {
        ds = IndexIJK(0, 1, 1);
        if (mr.Start.I == 0)
            d.I = 1;
    }
    else if (ms.J == 1)
    {
        ds = IndexIJK(1, 0, 1);
        if (mr.Start.J == 0)
            d.J = 1;
    }
    else if (ms.K == 1)
    {
        ds = IndexIJK(1, 1, 0);
        if (mr.Start.K == 0)
            d.K = 1;
    }
    else
    {
        assert(false);
    }

    cr.Start += d;
    cr.End += d;
    cr.Start += ds;

    return cr;
}

