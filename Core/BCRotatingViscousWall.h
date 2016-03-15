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
// $Id: BCViscousWall.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_BC_ROTATING_VISCOUS_WALL_H__
#define INCLUDED_BC_ROTATING_VISCOUS_WALL_H__

#include "BCPlanarLocal.h"
#include "BCViscousWall.h"
#include "RigidBodyMotion.h"

class BCRotatingViscousWall : public BCPlanarLocal
{
public:
    BCRotatingViscousWall(const IndexRange& meshRange, Direction direction, double TWall, const Vector3& angularVelocity);
    virtual ~BCRotatingViscousWall();

    virtual void SetMask(Structured<int>& mask) const;

    void ApplyTurMuK(Block& block)
    {
        BlockLocalFunctor ftor(block);
        Apply<BlockLocalFunctor>(block, ftor);
    }

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
    double mTWall;
    RotationalMotion* mRBM;
};

#endif // INCLUDED_BC_ROTATING_VISCOUS_WALL_H__

