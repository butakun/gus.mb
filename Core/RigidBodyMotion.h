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
// $Id: RigidBodyMotion.h 250 2012-07-20 01:52:56Z kato $
#ifndef INCLUDED_RIGID_BODY_MOTION_H__
#define INCLUDED_RIGID_BODY_MOTION_H__

#include "Vector3.h"
#include "Clock.h"
#include "Physics.h"
#include <cassert>

class RigidBodyMotion
{
public:
    virtual ~RigidBodyMotion() {}

    virtual RigidBodyMotion* Clone() const = 0;

    virtual Vector3 GetPosition(const Vector3& initial) const = 0;

protected:

private:
};

class TranslationalMotion : public RigidBodyMotion
{
public:
    TranslationalMotion(const Vector3& velocity) : mVelocity(velocity) {}
    TranslationalMotion(const TranslationalMotion& trans) : mVelocity(trans.mVelocity) {}
    virtual ~TranslationalMotion() {}

    TranslationalMotion* Clone() const { return new TranslationalMotion(*this); }

    virtual Vector3 GetPosition(const Vector3& initial) const
    {
        double t = Clock::GetInstance()->Time();
        return initial + t * mVelocity;
    }

protected:

private:
    Vector3 mVelocity;
};

class RotationalMotion : public RigidBodyMotion
{
public:
    RotationalMotion(const Vector3& origin, const Vector3& angularVelocity)
    : mOrigin(origin), mAngularVelocity(angularVelocity)
    {
        Initialize();
    }

    RotationalMotion(const RotationalMotion& rot)
    : mOrigin(rot.mOrigin), mAngularVelocity(rot.mAngularVelocity)
    {
        Initialize();
    }

    virtual ~RotationalMotion() {}

    virtual RotationalMotion* Clone() const { return new RotationalMotion(*this); }

    virtual Vector3 GetPosition(const Vector3& initial) const
    {
        double t = Clock::GetInstance()->Time();
        Vector3 p = initial - mOrigin;
        double a = dot_product(p, mAxis);
        Vector3 c = a * mAxis;
        Vector3 r = p - c;
        double rmag = r.Mag();
        Vector3 vt = cross_product(mAxis, p);
        Vector3 vr = cross_product(vt, mAxis);
        vt.Normalize();
        vr.Normalize();
        double theta = mOmega * t;
        double cos = std::cos(theta), sin = std::sin(theta);
        return mOrigin + c + rmag * cos * vr + rmag * sin * vt;
    }

    const Vector3& AngularVelocity() const { return mAngularVelocity; }
    const Vector3& Origin() const { return mOrigin; }
    double AngularSpeed() const { return mOmega; }
    const Vector3& Axis() const { return mAxis; }

    Vector3 RadialVector(const Vector3& pC) const
    {
        Vector3 vOC = pC - mOrigin;
        Vector3 vR = vOC - dot_product(vOC, mAxis) * mAxis;
        return vR;
    }

    Vector3 EntrainmentVelocityAt(const Vector3& pC) const
    {
        Vector3 vOC = pC - mOrigin;
        Vector3 v = cross_product(mAngularVelocity, vOC) / mVRef;
        return v;
    }

    Vector3 VelocityAt(const Vector3& pC) const
    {
        return EntrainmentVelocityAt(pC);
    }

protected:
    void Initialize()
    {
        mAxis = mAngularVelocity.Normalized();
        mOmega = mAngularVelocity.Mag();
        assert(mOmega > 0.0);
        mVRef = Physics::GetInstance()->VRef();
    }

private:
    Vector3 mOrigin;
    Vector3 mAngularVelocity;
    Vector3 mAxis;
    double mOmega;
    double mVRef; // cached copy of Physics::VRef
};

#endif // INCLUDED_RIGID_BODY_MOTION_H__

