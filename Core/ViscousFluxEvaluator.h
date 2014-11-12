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

