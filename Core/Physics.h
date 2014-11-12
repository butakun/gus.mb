// $Id: Physics.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_PHYSICS_H__
#define INCLUDED_PHYSICS_H__

#include "Sutherland.h"

/*
Nondimensionalization:

rho' = rho / rho0
u' = u / c0
v' = v / c0
w' = w / c0
p' = p / (rho0 * c0 * c0)
e' = e / (c0 * c0)
T' = T / (gamma * T0)
mu' = mu / mu0

t' = t * c0 / L

thus the equation of state becomes
p = rho * T

*/

class Physics
{
public:
    static void Initialize(
        double gamma = 1.4,
        double rhoRef = 1.225,
        double TRef = 288.15,
        double RGAS = 287.058
        );

    virtual ~Physics() {}

    double RGAS() const { return mRGAS; }
    double Gamma() const { return mGamma; }
    double PrandtlNumber() const { return mPr; }

    double LRef() const { return mLRef; }
    double TimeRef() const { return mLRef / mVRef; }

    double RhoRef() const { return mRhoRef; }
    double TRef() const { return mTRef; } // gamma * TRef
    double PRef() const { return mPRef; }
    double VRef() const { return mVRef; }
    double ERef() const { return mERef; } // specific total energy, unit = velocity squared
    double MuRef() const { return mMuRef; }
    double ReynoldsNumber() const { return mRhoRef * mVRef * mLRef / mMuRef; }

    double Mu(double T) const { return mViscosityModel.ViscosityNondimensional(T); }

    double TKERef() const { return mVRef * mVRef; }
    double OmegaRef() const { return mRhoRef * mVRef * mVRef / mMuRef; }

    const Sutherland& ViscosityModel() const { return mViscosityModel; }

    static Physics* GetInstance();

protected:
    Physics(double gamma, double rhoRef, double TRef, double RGAS);

private:
    static Physics* mPhysics;

    double mRGAS;
    double mGamma;
    double mPr;
    double mLRef;
    double mRhoRef, mVRef, mTRef, mPRef, mERef, mMuRef;

    Sutherland mViscosityModel;
};

#endif // INCLUDED_PHYSICS_H__

