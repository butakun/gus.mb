// $Id: IterationContext.h 262 2013-01-28 08:44:06Z kato $
#ifndef INCLUDED_ITERATION_CONTEXT_H__
#define INCLUDED_ITERATION_CONTEXT_H__

#include <iostream>

class IterationContext
{
public:
    enum IterationType { STEADY = 0, UNSTEADY = 1 };

    IterationContext(IterationType type = STEADY, double dt = -1.0);

    void SetType(IterationType type) { mIterationType = type; }

    bool IsUnsteady() const { return mIterationType == UNSTEADY; }

    void AdvanceTimeStep();
    void AdvanceInnerStep() { mInnerStep++; }

    int TimeStep() const { return mTimeStep; }
    int InnerStep() const { return mInnerStep; }
    double Time() const { return mT; }
    double DT() const { return mDT; }
    void SetDT(double dt) { mDT = dt; }

    double DTNonDim() const;

protected:

private:
    IterationType mIterationType;
    double mT;
    double mDT;
    int mTimeStep;
    int mInnerStep;
};

inline
std::ostream& operator << (std::ostream& o, const IterationContext& ic)
{
    o << "TimeStep: T = " << ic.Time() << ", " << ic.TimeStep() << "(" << ic.InnerStep() << ")";
    return o;
}

#endif // INCLUDED_ITERATION_CONTEXT_H__

