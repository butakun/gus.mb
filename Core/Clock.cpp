// $Id: Clock.cpp 41 2010-07-06 10:48:07Z kato $

#include "Clock.h"
#include <cstdlib>
#include <cassert>

Clock* Clock::mInstance = NULL;

void
Clock::Initialize(double t0)
{
    assert(Clock::mInstance == NULL);
    Clock::mInstance = new Clock(t0);
}

Clock*
Clock::GetInstance()
{
    if (Clock::mInstance == NULL)
    {
        Clock::mInstance = new Clock();
    }

    return Clock::mInstance;
}

Clock::Clock(double t0)
:   mT(t0)
{
}

double
Clock::Time() const
{
    return mT;
}

void
Clock::Advance(double dt)
{
    mT += dt;
}

