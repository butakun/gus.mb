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

