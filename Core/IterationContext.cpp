// $Id: IterationContext.cpp 262 2013-01-28 08:44:06Z kato $

#include "IterationContext.h"
#include "Physics.h"
#include <cassert>

IterationContext::IterationContext(IterationType type, double dt)
:   mIterationType(type), mT(0.0), mDT(dt), mTimeStep(0), mInnerStep(0)
{
    if (mIterationType == STEADY)
    {
        mDT = 0.0;
    }
}

void
IterationContext::AdvanceTimeStep()
{
    if (mIterationType == UNSTEADY)
    {
        assert(mDT > 0.0);
    }
    mT += mDT;
    mTimeStep++;
    mInnerStep = 0;
}

double
IterationContext::DTNonDim() const
{
    return mDT / Physics::GetInstance()->TimeRef();
}

