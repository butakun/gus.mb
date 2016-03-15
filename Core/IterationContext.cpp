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

