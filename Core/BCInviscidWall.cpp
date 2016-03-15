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
// $Id: BCInviscidWall.cpp 10 2010-04-04 03:12:58Z kato $

#if 0
#include "BCInviscidWall.h"
#include "Vector3.h"

BCInviscidWall::BCInviscidWall(const IndexRange& meshRange, Direction direction)
:   BCPlanarLocal(meshRange, direction)
{
}

BCInviscidWall::~BCInviscidWall()
{
}

void
BCInviscidWall::Apply(const Block& block, Structured<double>& U)
{
    IndexRange faceRange = MetricRange();
    Structured<double> surface = Surface(block);

    IndexIJK dGhost(0, 0, 0), dInterior(0, 0, 0);
    switch (PatchDirection())
    {
    case I:
        dInterior.I = 1;
        break;
    case INEG:
        dGhost.I = 1;
        break;
    case J:
        dInterior.J = 1;
        break;
    case JNEG:
        dGhost.J = 1;
        break;
    case K:
        dInterior.K = 1;
        break;
    case KNEG:
        dGhost.K = 1;
        break;
    default:
        assert(false);
    }

    for (int k = faceRange.Start.K; k <= faceRange.End.K; ++k)
    {
        for (int j = faceRange.Start.J; j <= faceRange.End.J; ++j)
        {
            for (int i = faceRange.Start.I; i <= faceRange.End.I; ++i)
            {
                IndexIJK iFace(i, j, k);
                IndexIJK iGhost = iFace + dGhost;
                IndexIJK iInterior = iFace + dInterior;

                Vector3 Sn(surface(iFace));

                double* UGhost = U(iGhost);
                double* UInterior = U(iInterior);

                for (int l = 0; l < 5; ++l)
                {
                    UGhost[l] = UInterior[l];
                }

                Vector3 V(UInterior[1], UInterior[2], UInterior[3]);
                V = V / UInterior[0];

                Vector3 n = Sn / Sn.Mag();
                Vector3 Vn = n * dot_product(V, n);
                Vector3 Vghost = V - 2.0 * Vn;
                UGhost[1] = UGhost[0] * Vghost[0];
                UGhost[2] = UGhost[1] * Vghost[1];
                UGhost[3] = UGhost[2] * Vghost[2];
            }
        }
    }
}
#endif

