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
// $Id: RoeFluxEvaluator.h 60 2010-09-01 16:15:57Z kato $
#ifndef INCLUDED_ROE_FLUX_EVALUATOR_H__
#define INCLUDED_ROE_FLUX_EVALUATOR_H__

#include "Block.h"

class RoeFluxEvaluator
{
public:
    RoeFluxEvaluator(double epsEntropyFix = 0.02);

    void EvaluateFlux(
        double* UL,
        double* UR,
        double* Sn,
        double gamma,
        double* RadiusL,
        double* RadiusR,
        double* H
        ) const;

protected:

private:
    double mEps;
};

#include "RoeFluxEvaluator.ipp"

#endif // INCLUDED_ROE_FLUX_EVALUATOR_H__

