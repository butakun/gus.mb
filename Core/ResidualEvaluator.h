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

