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
#ifndef INCLUDED_BC_PLANAR_H__
#define INCLUDED_BC_PLANAR_H__

#include "BC.h"
#include "IndexUtils.h"

class BCPlanar : public BC
{
public:
    BCPlanar(const IndexRange& meshRange, Direction direction);
    virtual ~BCPlanar() {}

    Direction PatchDirection() const { return mDirection; }
    IndexRange MetricRange() const { return IndexUtils::FromMeshRangeToMetricRange(MeshRange(), PatchDirection()); }
    IndexRange RindRange() const;
    Structured<double> Surface(const Block& block) const
    {
        switch (PatchDirection())
        {
        case I:
        case INEG:
            return block.Sxi();
        case J:
        case JNEG:
            return block.Seta();
        case K:
        case KNEG:
            return block.Szeta();
        }
        assert(false);
        return block.Sxi();
    }

protected:
    IndexIJK DGhost() const { return mDGhost; }
    IndexIJK DInterior() const { return mDInterior; }
    IndexIJK DeltaInterior() const { return mDeltaInterior; }
    double SnDirection() const { return mSnDirection; }
    void SetMaskWithAValue(Structured<int>& mask, int value) const;

private:
    Direction mDirection;
    IndexIJK mDGhost, mDInterior, mDeltaInterior;
    double mSnDirection; // 1.0 if the face metric (Sn) is pointed inward (from boundary toward interior cells), -1.0 if outward.
};

#endif // INCLUDED_BC_PLANAR_H__

