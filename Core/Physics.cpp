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
// $Id: Physics.cpp 255 2012-12-05 07:47:00Z kato $

#include "Physics.h"
#include <cstddef>
#include <cmath>
#include <cassert>

Physics* Physics::mPhysics = NULL;

Physics*
Physics::GetInstance()
{
    assert(Physics::mPhysics != NULL);
    return Physics::mPhysics;
}

void
Physics::Initialize(double gamma, double rhoRef, double TRef, double RGAS)
{
    assert(mPhysics == NULL);
    Physics::mPhysics = new Physics(gamma, rhoRef, TRef, RGAS);
}

Physics::Physics(double gamma, double rhoRef, double TRef, double RGAS)
:   mRGAS(RGAS), mGamma(gamma), mPr(0.72), mRhoRef(rhoRef), mTRef(mGamma * TRef),
    mViscosityModel(120.0, 291.15, 18.27e-6, mTRef)
{
    mLRef = 1.0;
    mERef = mGamma * mRGAS * TRef;
    mVRef = std::sqrt(mERef);
    mPRef = mRhoRef * mERef;
    mMuRef = mViscosityModel.ViscosityDimensional(TRef);
}

