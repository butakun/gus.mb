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
// $Id: BC.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_BC_H__
#define INCLUDED_BC_H__

#include "Block.h"
#include "Structured.h"

class BC
{
public:
    BC(const IndexRange& meshRange);
    virtual ~BC() {}

    virtual void Apply(const Block& block, Structured<double>& U) = 0;
    virtual void ApplyTurb(const Block& block, Structured<double>& UT, const Structured<double>& U) = 0;
    virtual void SetMask(Structured<int>& mask) const {}

    virtual void PreIteration() {}
    virtual void PostIteration() {}

    const IndexRange& MeshRange() const { return mMeshRange; }

protected:

private:
    IndexRange mMeshRange;
};

#endif // INCLUDED_BC_H__

