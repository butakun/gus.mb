// $Id: SolverUtils.h 74 2010-09-20 17:33:05Z kato $
#ifndef INCLUDED_SOLVER_UTILS_H__
#define INCLUDED_SOLVER_UTILS_H__

#include "Block.h"
#include "Structured.h"

namespace SolverUtils
{

void ComputeSpectralRadius(
    const Block& block,
    const Structured<double>& U,
    Structured<double>& Lambda
    );

double ComputeTimeStep(
    const Block& block,
    const Structured<double>& U,
    Structured<double>& DT,
    double cfl,
    bool localStepping = true
    );

}

#endif // INCLUDED_SOLVER_UTILS_H__

