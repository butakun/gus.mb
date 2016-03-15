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
// $Id: ViscousFluxEvaluator.h 52 2010-07-23 08:53:22Z kato $
#ifndef INCLDUED_VISCOUS_FLUX_EVALUATOR_H__
#define INCLDUED_VISCOUS_FLUX_EVALUATOR_H__

#include "Structured.h"
class Block;

class ViscousFluxEvaluator
{
public:
    ViscousFluxEvaluator();

    void EvaluateFlux(
        const Structured<double>& R,
        const Block& block,
        const Structured<double>& U
        ) const;

protected:

private:
};

#endif // INCLDUED_VISCOUS_FLUX_EVALUATOR_H__

