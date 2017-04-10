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

