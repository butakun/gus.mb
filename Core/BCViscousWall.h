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
#ifndef INCLUDED_BC_VISCOUS_WALL_H__
#define INCLUDED_BC_VISCOUS_WALL_H__

#include "BCPlanarLocal.h"

class BlockLocalFunctor
{
public:
    BlockLocalFunctor(Block& block) : mBlock(block) {}
    void operator () (
        const IndexIJK& iFace,
        const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
        const Vector3& Sn, const Block& block
        )
    {
        Structured<double> TurMuK = mBlock.TurMuK();
        double* turmukI = TurMuK(iInterior);
        double* turmukG = TurMuK(iGhost);
        turmukG[0] = -turmukI[0];
        turmukG[1] = -turmukI[1];
    }

private:
    Block& mBlock;
};

class BCViscousWall : public BCPlanarLocal
{
public:
    BCViscousWall(const IndexRange& meshRange, Direction direction, bool isStationary = false, double TWall = -1.0);
    virtual ~BCViscousWall() {}

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
    bool mIsStationary;
    double mTWall;
};

#endif // INCLUDED_BC_VISCOUS_WALL_H__

