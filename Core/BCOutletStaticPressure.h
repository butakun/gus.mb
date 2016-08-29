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
#ifndef INCLUDE_BC_OUTLET_STATIC_PRESSURE_H__
#define INCLUDE_BC_OUTLET_STATIC_PRESSURE_H__

#include "BCPlanarLocal.h"

class BCOutletStaticPressure : public BCPlanarLocal
{
public:
    BCOutletStaticPressure(const IndexRange& meshRange, Direction direction, double pressure, bool forcePressure = false);
    virtual ~BCOutletStaticPressure() {}

protected:
    virtual void LocalFunc(
        const IndexIJK& iFace,
        const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
        const Vector3& Sn, Structured<double>& U, const Block& block
        );

    virtual void LocalFuncTurb(
        const IndexIJK& iFace,
        const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
        const Vector3& Sn, Structured<double>& UT, const Structured<double>& U, const Block& block
        );

private:
    double mStaticPressure;
    bool mForcePressure;
};

#endif // INCLUDE_BC_OUTLET_STATIC_PRESSURE_H__

