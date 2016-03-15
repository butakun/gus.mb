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

