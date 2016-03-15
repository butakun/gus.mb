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

