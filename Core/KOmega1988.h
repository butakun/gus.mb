// $Id: KOmega1988.h 284 2013-06-14 03:22:13Z kato $
#ifndef INCLUDED_K_OMEGA_1988_H__
#define INCLUDED_K_OMEGA_1988_H__

#include "Model.h"
#include "Vector3.h"
#include "Structured.h"
#include "Block.h"
#include "TurbulenceSpec.h"

class Residual;
class AuxiliaryCell;

namespace KOmega1988
{

const double SIGMA_K = 0.5, SIGMA_OMEGA = 0.5, BETA_K = 0.09, BETA_OMEGA = 0.075, ALPHA_OMEGA = 5.0 / 9.0;

class ConvectiveFluxEvaluator
{
public:
    ConvectiveFluxEvaluator() {}
    virtual ~ConvectiveFluxEvaluator() {}

    virtual void EvaluateFlux(double* H, double* UL, double* UR, double* UTL, double* UTR, const Vector3& Sn);

protected:

private:
};

class DiffusiveFluxEvaluator
{
public:
    DiffusiveFluxEvaluator(double Re) : mRe(Re) {}
    virtual ~DiffusiveFluxEvaluator() {}

    virtual void EvaluateFlux(
        double* H,
        const AuxiliaryCell& cell,
        const Structured<double>& UT,
        const Structured<double>& U,
        const Structured<double>& MuK,
        const Structured<double>& TurMuK,
        const Vector3& Sn
        );

protected:

private:
    double mRe;
};

class ResidualEvaluator
{
public:
    ResidualEvaluator();
    virtual ~ResidualEvaluator() {}

    Residual EvaluateResidual(
        const Block& block,
        const Structured<double>& UT,
        Structured<double>& R
        );

protected:

private:
};

class KOmegaAdaptor
{
public:
    KOmegaAdaptor(const Structured<double>& U_, const Structured<double>& UT_) : U(U_), UT(UT_) {}
    void Evaluate(double* UU, const IndexIJK& i) const
    {
        UU[0] = UT(i)[0] / U(i)[0];
        UU[1] = UT(i)[1] / U(i)[0];
    }

protected:

private:
    const Structured<double>& U;
    const Structured<double>& UT;
};

class TurbulenceModel : public Model
{
public:
    TurbulenceModel();

    int DOF() const { return 2; }

    double MuT(const double* U, const double* UT) const;

    Structured<double>& GetU(Block& block) const { return block.UT(); }
    const Structured<double>& GetU(const Block& block) const { return block.UT(); }

    void ApplyBCs(Block& block) const;

    double ComputeTimeStep(
        const Block& block,
        const Structured<double>& U,
        Structured<double>& DT,
        double cfl,
        bool localTimeStepping
        );

    // For matrix-free LUSGS
    double ScalarCoeff(
        const Block& block, const Structured<double>& U,
        const IndexIJK& Ii, const IndexIJK& Ij,
        double* Sij, double SijSign, double SijAbs, const Structured<double>& Radius
        ) const;
    void JacobianDU(
        double* JacDU, const Block& block, const IndexIJK& i,
        double* U, double* dU, double* Sij, double SijSign, double SijAbs, double* Radius
        ) const;
#define SCALAR_DIAG 0
#if SCALAR_DIAG
    double DiagonalAddition(const Block& block, const Structured<double>& UT, const IndexIJK& Ii) const;
#else
    void DiagonalAddition(double* D, const Block& block, const Structured<double>& UT, const IndexIJK& Ii) const;
#endif

    // Coordinate frame transform
    virtual void FromGlobalToLocal(double* ULocal, double* UGlobal, const Block& block, const IndexIJK& ijk) const
    {
        ULocal[0] = UGlobal[0]; ULocal[1] = UGlobal[1];
    }

    virtual void FromLocalToGlobal(double* UGlobal, double* ULocal, const Block& block, const IndexIJK& ijk) const
    {
        UGlobal[0] = ULocal[0]; UGlobal[1] = ULocal[1];
    }

    void EvaluateUT(double* UT, const double* U, const TurbulenceSpec& turbSpec) const
    {
        double tke = turbSpec.TKE_Nondimensional();
        double omega = turbSpec.Omega_Nondimensional();
        double rho = U[0];
        UT[0] = rho * tke;
        UT[1] = rho * omega;
    }

    void SetUT(Structured<double>& UT, const Structured<double>& U, const IndexRange& range, const TurbulenceSpec& turbSpec) const;

protected:

private:
};

}

#include "KOmega1988.ipp"

#endif // INCLUDED_K_OMEGA_1988_H__

