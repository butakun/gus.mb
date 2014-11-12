// $Id: TurbulenceSpec.cpp 280 2013-06-06 09:24:16Z kato $

#include "TurbulenceSpec.h"
#include "Physics.h"
#include <cassert>

std::ostream& operator << (std::ostream& o, const TurbulenceSpec& spec)
{
    spec.Dump(o);
    return o;
}

double
TurbulenceSpec::TKE_Nondimensional() const
{
    double VRef = Physics::GetInstance()->VRef();

    return TKE_Dimensional() / (VRef * VRef);
}

double
TurbulenceSpec::Omega_Nondimensional() const
{
    Physics* phys = Physics::GetInstance();
    double muRef = phys->MuRef();
    double rhoRef = phys->RhoRef();
    double velRef = phys->VRef();

    return Omega_Dimensional() * muRef / (rhoRef * velRef * velRef);
}

TurbulenceSpecIntensityViscosityRatio::TurbulenceSpecIntensityViscosityRatio(
    double I,
    double viscRatio,
    double rho, // reference density (dimensional), used for evaluation of omega from viscRatio
    double U, // reference velocity (dimensional) on which turbulence intensity is based.
    double T  // reference temperature (dimensional) on which viscosity is based.
    )
:   mI(I), mViscRatio(viscRatio), mRho(rho), mU(U), mT(T)
{
    assert(mRho > 0.0 && mU > 0.0 && mT > 0.0);
}

TurbulenceSpecIntensityViscosityRatio*
TurbulenceSpecIntensityViscosityRatio::Clone() const
{
    return new TurbulenceSpecIntensityViscosityRatio(*this);
}

void
TurbulenceSpecIntensityViscosityRatio::Dump(std::ostream& o) const
{
    o << "IntensityViscosityRatio: I = " << mI << ", ViscRatio = " << mViscRatio << ", Rho = " << mRho << ", U = " << mU << ", T = " << mT;
}

double
TurbulenceSpecIntensityViscosityRatio::TKE_Dimensional() const
{
    return 1.5 * (mI * mU) * (mI * mU);
}

double
TurbulenceSpecIntensityViscosityRatio::Omega_Dimensional() const
{
    double mu = Physics::GetInstance()->ViscosityModel().ViscosityDimensional(mT);
    double mut = mViscRatio * mu;
    return mRho * TKE_Dimensional() / mut;
}

TurbulenceSpecKOmegaNondimensional::TurbulenceSpecKOmegaNondimensional(double k_nondim, double omega_nondim)
:   mTKENondim(k_nondim), mOmegaNondim(omega_nondim)
{
}

double
TurbulenceSpecKOmegaNondimensional::TKE_Dimensional() const
{
    Physics* phys = Physics::GetInstance();
    double vref = phys->VRef();
    return mTKENondim * vref * vref;
}

double
TurbulenceSpecKOmegaNondimensional::TKE_Nondimensional() const
{
    return mTKENondim;
}

double
TurbulenceSpecKOmegaNondimensional::Omega_Dimensional() const
{
    Physics* phys = Physics::GetInstance();
    double rhoref = phys->RhoRef();
    double muref = phys->MuRef();
    double vref = phys->VRef();
    return mOmegaNondim * rhoref * vref * vref / muref;
}

double
TurbulenceSpecKOmegaNondimensional::Omega_Nondimensional() const
{
    return mOmegaNondim;
}

