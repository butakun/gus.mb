// $Id: TurbulenceSpec.h 280 2013-06-06 09:24:16Z kato $
#ifndef INCLUDED_TURBULENCE_SPEC_H__
#define INCLUDED_TURBULENCE_SPEC_H__

#include <iostream>

class TurbulenceSpec
// represents a turbulence state independent of turbulence models
{
public:
    virtual ~TurbulenceSpec() {}

    virtual TurbulenceSpec* Clone() const = 0;

    virtual void Dump(std::ostream& o) const = 0;

    virtual double TKE_Dimensional() const = 0;
    virtual double Omega_Dimensional() const = 0;

    virtual double TKE_Nondimensional() const;
    virtual double Omega_Nondimensional() const;

protected:

private:
};

std::ostream& operator << (std::ostream& o, const TurbulenceSpec& spec);

class NoTurbulence : public TurbulenceSpec
{
public:
    NoTurbulence() {}
    virtual ~NoTurbulence() {}

    virtual NoTurbulence* Clone() const { return new NoTurbulence(*this); }

    virtual void Dump(std::ostream& o) const { o << "NoTurbulence"; }

    virtual double TKE_Dimensional() const { return 0.0; }
    virtual double Omega_Dimensional() const { return 1.0; }

protected:

private:
};

class TurbulenceSpecIntensityViscosityRatio : public TurbulenceSpec
{
public:
    TurbulenceSpecIntensityViscosityRatio(
        double I,
        double viscRatio,
        double rho, // reference density (dimensional), used for evaluation of omega from viscRatio
        double U, // reference velocity (dimensional) on which turbulence intensity is based.
        double T // reference temperature (dimensional) on which viscosity is based.
        );

    virtual ~TurbulenceSpecIntensityViscosityRatio() {}

    virtual TurbulenceSpecIntensityViscosityRatio* Clone() const;

    virtual void Dump(std::ostream& o) const;

    virtual double TKE_Dimensional() const;
    virtual double Omega_Dimensional() const;

protected:

private:
    double mI, mViscRatio, mRho, mU, mT;
};

class TurbulenceSpecKOmegaNondimensional : public TurbulenceSpec
{
public:
    TurbulenceSpecKOmegaNondimensional(double k_nondim, double omega_nondim);
    virtual ~TurbulenceSpecKOmegaNondimensional() {}

    virtual TurbulenceSpecKOmegaNondimensional* Clone() const { return new TurbulenceSpecKOmegaNondimensional(*this); }

    virtual void Dump(std::ostream& o) const { o << "KOmegaNondimensional: tke = " << mTKENondim << ", omega = " << mOmegaNondim; }

    virtual double TKE_Dimensional() const;
    virtual double TKE_Nondimensional() const;
    virtual double Omega_Dimensional() const;
    virtual double Omega_Nondimensional() const;

protected:

private:
    double mTKENondim, mOmegaNondim;
};

#endif // INCLUDED_TURBULENCE_SPEC_H__

