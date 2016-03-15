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
// $Id: BCEnforce2D.h 175 2012-01-04 06:01:45Z kato $
#ifndef INCLUDED_ENFORCE_2D_H__
#define INCLUDED_ENFORCE_2D_H__

#include "BC.h"

class BCEnforce2D : public BC
{
public:
    enum Type { TWOD_X, TWOD_Y, TWOD_Z, AXISYM_X, AXISYM_Y, AXISYM_Z };

    BCEnforce2D(Type type);
    virtual ~BCEnforce2D();

    virtual void Apply(const Block& block, Structured<double>& U);
    virtual void ApplyTurb(const Block& block, Structured<double>& UT);

protected:

private:
    Type mType;
};

#endif // INCLUDED_ENFORCE_2D_H__

