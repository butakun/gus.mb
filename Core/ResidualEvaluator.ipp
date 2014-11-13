// $Id: ResidualEvaluator.ipp 258 2012-12-26 07:36:54Z kato $

#include "Physics.h"
#include "SolverFunctions.h"
#include "Reconstructor.h"
#include "RoeFluxEvaluator.h"
#include "ViscousFluxEvaluator.h"
#include "Residual.h"
#include "RigidBodyMotion.h"
#include <cstdlib>

inline void ViscousWallFlux(double* UL, double* UR, double* Sn, double gamma, double* H)
{
    double U[5];
    for (int l = 0; l < 5; ++l)
        U[l] = 0.5 * (UL[l] + UR[l]);
    double rhoe, p;
    rhoe = U[4] - 0.5 * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) / U[0]; // FIXME: rothalpy
    p = (gamma - 1.0) * rhoe;
    H[0] = 0.0;
    H[1] = Sn[0] * p;
    H[2] = Sn[1] * p;
    H[3] = Sn[2] * p;
    H[4] = 0.0;
}

template<int N, class RECON>
inline
ResidualEvaluator<N, RECON>::ResidualEvaluator()
{
}

template<int N, class RECON>
inline
void
ResidualEvaluator<N, RECON>::EvaluateResidual(
    const Block& block,
    const Structured<double>& U,
    Structured<double>& R
    ) const
{
    IndexRange cr = block.CellRange();

    const Structured<double>& Sxi = block.Sxi();
    const Structured<double>& Seta = block.Seta();
    const Structured<double>& Szeta = block.Szeta();
    const Structured<double>& Vol = block.Vol();
    const Structured<double>& Radius = block.Radius();
    const Structured<int>& Mask = block.Mask();

    double gamma = Physics::GetInstance()->Gamma();
    double H[N];
    R.SetTo(0.0);

    RoeFluxEvaluator fluxEval;
    ViscousFluxEvaluator viscEval;

    bool ViscWall = true;

#define MUSCL 1

    // Xi flux
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I - 1; i <= cr.End.I; ++i)
            {
                double UL[N], UR[N], *Sn, *RadL, *RadR;
                double *RL, *RR;
#if MUSCL
                mRecon.Reconstruct(UL, UR, U, Radius, IndexIJK(i, j, k), IndexIJK(1, 0, 0));
#else
                for (int l = 0; l < N; ++l)
                {
                    UL[l] = U(i, j, k)[l];
                    UR[l] = U(i + 1, j, k)[l];
                }
#endif
                Sn = Sxi(i, j, k);
                RadL = Radius(i, j, k);
                RadR = Radius(i + 1, j, k);
                if (ViscWall && (*Mask(i, j, k) == 1 || *Mask(i + 1, j, k) == 1))
                {
                    ViscousWallFlux(UL, UR, Sn, gamma, H);
                }
                else
                {
                    fluxEval.EvaluateFlux(UL, UR, Sn, gamma, RadL, RadR, H);
                }

                if (i > cr.Start.I - 1)
                {
                    RL = R(i, j, k);
                    for (int l = 0; l < N; ++l)
                    {
                        RL[l] += H[l];
                    }
                }
                if (i < cr.End.I)
                {
                    RR = R(i + 1, j, k);
                    for (int l = 0; l < N; ++l)
                    {
                        RR[l] -= H[l];
                    }
                }
            }
        }
    }

    // Eta flux
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J - 1; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double UL[N], UR[N], *Sn, *RadL, *RadR;
                double *RL, *RR;
#if MUSCL
                mRecon.Reconstruct(UL, UR, U, Radius, IndexIJK(i, j, k), IndexIJK(0, 1, 0));
#else
                for (int l = 0; l < N; ++l)
                {
                    UL[l] = U(i, j, k)[l];
                    UR[l] = U(i, j + 1, k)[l];
                }
#endif
                Sn = Seta(i, j, k);
                RadL = Radius(i, j, k);
                RadR = Radius(i, j + 1, k);
                if (ViscWall && (*Mask(i, j, k) == 1 || *Mask(i, j + 1, k) == 1))
                {
                    ViscousWallFlux(UL, UR, Sn, gamma, H);
                }
                else
                {
                    fluxEval.EvaluateFlux(UL, UR, Sn, gamma, RadL, RadR, H);
                }

                if (j > cr.Start.J - 1)
                {
                    RL = R(i, j, k);
                    for (int l = 0; l < N; ++l)
                    {
                        RL[l] += H[l];
                    }
                }
                if (j < cr.End.J)
                {
                    RR = R(i, j + 1, k);
                    for (int l = 0; l < N; ++l)
                    {
                        RR[l] -= H[l];
                    }
                }
            }
        }
    }

    // Zeta flux
    for (int k = cr.Start.K - 1; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double UL[N], UR[N], *Sn, *RadL, *RadR;
                double *RL, *RR;
#if MUSCL
                mRecon.Reconstruct(UL, UR, U, Radius, IndexIJK(i, j, k), IndexIJK(0, 0, 1));
#else
                for (int l = 0; l < N; ++l)
                {
                    UL[l] = U(i, j, k)[l];
                    UR[l] = U(i, j, k + 1)[l];
                }
#endif
                Sn = Szeta(i, j, k);
                RadL = Radius(i, j, k);
                RadR = Radius(i, j, k + 1);
                if (ViscWall && (*Mask(i, j, k) == 1 || *Mask(i, j, k + 1) == 1))
                {
                    ViscousWallFlux(UL, UR, Sn, gamma, H);
                }
                else
                {
                    fluxEval.EvaluateFlux(UL, UR, Sn, gamma, RadL, RadR, H);
                }

                if (k > cr.Start.K - 1)
                {
                    RL = R(i, j, k);
                    for (int l = 0; l < N; ++l)
                    {
                        RL[l] += H[l];
                    }
                }
                if (k < cr.End.K)
                {
                    RR = R(i, j, k + 1);
                    for (int l = 0; l < N; ++l)
                    {
                        RR[l] -= H[l];
                    }
                }
            }
        }
    }

    // Viscous Flux
    viscEval.EvaluateFlux(R, block, U);

    // Source Terms
    if (block.IsRotating())
    {
        RotationalMotion* rm = dynamic_cast<RotationalMotion*>(block.GetRigidBodyMotion());
        assert(rm != NULL);
        Vector3 omega = rm->AngularVelocity();

        const Physics* PHYS = Physics::GetInstance();
        double Fref = PHYS->VRef() * PHYS->VRef() / PHYS->LRef();

        for (int k = cr.Start.K; k <= cr.End.K; ++k)
        {
            for (int j = cr.Start.J; j <= cr.End.J; ++j)
            {
                for (int i = cr.Start.I; i <= cr.End.I; ++i)
                {
                    double* rhs = R(i, j, k);
                    double vol = *Vol(i, j, k);
                    double* rad = Radius(i, j, k);
                    Vector3 r(rad[0], rad[1], rad[2]);
                    double* UU = U(i, j, k);
                    double rho = UU[0];
                    Vector3 v(UU[1] / rho, UU[2] / rho , UU[3] / rho);
                    Vector3 Fcentrifugal = -cross_product(omega, cross_product(omega, r));
                    Vector3 Fcoriolis = -2.0 * cross_product(omega, v);
                    rhs[1] += vol * rho / Fref * (Fcentrifugal.X() + Fcoriolis.X());
                    rhs[2] += vol * rho / Fref * (Fcentrifugal.Y() + Fcoriolis.Y());
                    rhs[3] += vol * rho / Fref * (Fcentrifugal.Z() + Fcoriolis.Z());
                }
            }
        }
    }
}

