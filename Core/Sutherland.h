// $Id: Sutherland.h 33 2010-06-01 17:29:09Z kato $
#ifndef INCLUDED_SUTHERLAND_H__
#define INCLUDED_SUTHERLAND_H__

#include <cmath>

class Sutherland
{
public:
    // TRef is the reference temperature of the problem,
    // not the temperature (T0) corresponding to the reference viscosity (Mu0)
    // in the definition of the Sutherland's formula.
    Sutherland(double C, double T0, double Mu0, double TRef)
    : mC(C), mT0(T0), mMu0(Mu0), mCNondim(mC / TRef)
    {}

    double ViscosityDimensional(double T) const
    {
        return mMu0 * (mT0 + mC) / (T + mC) * std::pow(T / mT0, 1.5);
    }

    double ViscosityNondimensional(double T) const
    {
        return (1.0 + mCNondim) / (T + mCNondim) * std::pow(T, 1.5);
    }

protected:

private:
    double mC, mT0, mMu0;
    double mCNondim;
};

#endif // INCLUDED_SUTHERLAND_H__

