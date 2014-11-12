// $Id: BCFunctors.cpp 4 2010-03-06 14:10:00Z kato $

#include "BCFunctors.h"
#include <cmath>

BCInviscidWall::BCInviscidWall(const IndexRange& range, int direction)
:   mRange(range), mDirection(direction)
{
}

BCInviscidWall::~BCInviscidWall()
{
}

void
BCInviscidWall::Apply(
    Structured<double> U,
    const Block& block
    )
{
    for (int i = mRange.Start.I; i <= mRange.End.I; ++i)
    {
        double* sn = block.Seta()(i, mJWall, 0);
        double* Ui = block.U()(i, mJIn, 0); // interior cell
        double* Ug = block.U()(i, mJGhost, 0); // ghost cell
        double snAbs, v[3], vn, vg[3];

        snAbs = std::sqrt(sn[0] * sn[0] + sn[1] * sn[1] + sn[2] * sn[2]);

        v[0] = Ui[1] / Ui[0];
        v[1] = Ui[2] / Ui[0];
        v[2] = Ui[3] / Ui[0];

        vn = (v[0] * sn[0] + v[1] * sn[1] + v[2] * sn[2]) / snAbs;

        vg[0] = v[0] - 2.0 * vn * sn[0] / snAbs;
        vg[1] = v[1] - 2.0 * vn * sn[1] / snAbs;
        vg[2] = v[2] - 2.0 * vn * sn[2] / snAbs;

        Ug[0] = Ui[0];
        Ug[1] = Ui[0] * vg[0];
        Ug[2] = Ui[0] * vg[1];
        Ug[3] = Ui[0] * vg[2];
        Ug[4] = Ui[4];
    }
}

BCExtrapolate::BCExtrapolate(int iin, int ighost)
:   mIIn(iin), mIGhost(ighost)
{
}

BCExtrapolate::~BCExtrapolate()
{
}

void
BCExtrapolate::Apply(
    Structured<double> U,
    const Block& block
    )
{
    for (int j = block.CellRange().Start.J; j < block.CellRange().End.J; ++j)
    {
        for (int k = block.CellRange().Start.K; k < block.CellRange().End.K; ++k)
        {
            double* Ui = U(mIIn, j, k);
            double* Ug = U(mIGhost, j, k);
            for (int l = 0; l < 5; ++l)
            {
                Ug[l] = Ui[l];
            }
        }
    }
}

