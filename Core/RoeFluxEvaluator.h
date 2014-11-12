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

