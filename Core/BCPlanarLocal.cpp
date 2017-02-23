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

#include "BCPlanarLocal.h"
#include "Vector3.h"

void
BCPlanarLocal::Apply(const Block& block, Structured<double>& U)
{
    IndexRange faceRange = MetricRange();
    Structured<double> surface = Surface(block);

    IndexIJK dGhost = DGhost();
    IndexIJK dInterior = DInterior();
    IndexIJK deltaInterior = DeltaInterior();
    double SnDir = SnDirection();
#if 0
    switch (PatchDirection())
    {
    case I:
        dInterior.I = 1;
        deltaInterior.I = 1;
        SnDir = 1.0;
        break;
    case INEG:
        dGhost.I = 1;
        deltaInterior.I = -1;
        SnDir = -1.0;
        break;
    case J:
        dInterior.J = 1;
        deltaInterior.J = 1;
        SnDir = 1.0;
        break;
    case JNEG:
        dGhost.J = 1;
        deltaInterior.J = -1;
        SnDir = -1.0;
        break;
    case K:
        dInterior.K = 1;
        deltaInterior.K = 1;
        SnDir = 1.0;
        break;
    case KNEG:
        dGhost.K = 1;
        deltaInterior.K = -1;
        SnDir = -1.0;
        break;
    default:
        assert(false);
    }
#endif

    for (int k = faceRange.Start.K; k <= faceRange.End.K; ++k)
    {
        for (int j = faceRange.Start.J; j <= faceRange.End.J; ++j)
        {
            for (int i = faceRange.Start.I; i <= faceRange.End.I; ++i)
            {
                IndexIJK iFace(i, j, k);
                IndexIJK iGhost = iFace + dGhost;
                IndexIJK iInterior = iFace + dInterior;

                Vector3 Sn = SnDir * Vector3(surface(iFace));

                LocalFunc(iFace, iGhost, iInterior, deltaInterior, Sn, U, block);
            }
        }
    }
}

void
BCPlanarLocal::ApplyTurb(const Block& block, Structured<double>& UT, const Structured<double>& U)
{
    IndexRange faceRange = MetricRange();
    Structured<double> surface = Surface(block);

    IndexIJK dGhost = DGhost();
    IndexIJK dInterior = DInterior();
    IndexIJK deltaInterior = DeltaInterior();
    double SnDir = SnDirection();
#if 0
    IndexIJK dGhost(0, 0, 0), dInterior(0, 0, 0), deltaInterior(0, 0, 0);
    switch (PatchDirection())
    {
    case I:
        dInterior.I = 1;
        deltaInterior.I = 1;
        break;
    case INEG:
        dGhost.I = 1;
        deltaInterior.I = -1;
        break;
    case J:
        dInterior.J = 1;
        deltaInterior.J = 1;
        break;
    case JNEG:
        dGhost.J = 1;
        deltaInterior.J = -1;
        break;
    case K:
        dInterior.K = 1;
        deltaInterior.K = 1;
        break;
    case KNEG:
        dGhost.K = 1;
        deltaInterior.K = -1;
        break;
    default:
        assert(false);
    }
#endif

    for (int k = faceRange.Start.K; k <= faceRange.End.K; ++k)
    {
        for (int j = faceRange.Start.J; j <= faceRange.End.J; ++j)
        {
            for (int i = faceRange.Start.I; i <= faceRange.End.I; ++i)
            {
                IndexIJK iFace(i, j, k);
                IndexIJK iGhost = iFace + dGhost;
                IndexIJK iInterior = iFace + dInterior;

                Vector3 Sn = SnDir * Vector3(surface(iFace));

                LocalFuncTurb(iFace, iGhost, iInterior, deltaInterior, Sn, UT, U, block);
            }
        }
    }
}

