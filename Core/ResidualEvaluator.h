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
// $Id: ResidualEvaluator.h 244 2012-06-01 15:39:32Z kato $
#ifndef INCLUDED_RESIDUAL_EVALUATOR_H__
#define INCLUDED_RESIDUAL_EVALUATOR_H__

#include "Communicator.h"
#include "Block.h"

class ResidualEvaluatorBase
{
public:
    virtual ~ResidualEvaluatorBase() {}

    virtual void EvaluateResidual(
        const Block& block,
        const Structured<double>& U,
        Structured<double>& R
        ) const = 0;
};

template <int N, class RECON>
class ResidualEvaluator : public ResidualEvaluatorBase
{
public:
    ResidualEvaluator();
    virtual ~ResidualEvaluator() {}

    virtual void EvaluateResidual(
        const Block& block,
        const Structured<double>& U,
        Structured<double>& R
        ) const;

protected:

private:
    RECON mRecon;
};

#include "ResidualEvaluator.ipp"

#endif // INCLUDED_RESIDUAL_EVALUATOR_H__

