// $Id: KOmega1988.cpp 284 2013-06-14 03:22:13Z kato $

#include "KOmega1988.h"
#include "AuxiliaryCell.h"
#include "Residual.h"
#include "GradientEvaluator.h"
#include <cassert>

#define DEBUG 0

class VelocityAdaptor
{
public:
    VelocityAdaptor(const Structured<double>& UU_) : UU(UU_) {}

    void Evaluate(double* U, const IndexIJK& i) const
    {
        double* uu = UU(i);
        U[0] = uu[1] / uu[0];
        U[1] = uu[2] / uu[0];
        U[2] = uu[3] / uu[0];
    }

private:
    const Structured<double>& UU;
};

using namespace KOmega1988;

ResidualEvaluator::ResidualEvaluator()
{
}

Residual
ResidualEvaluator::EvaluateResidual(
    const Block& block,
    const Structured<double>& UT,
    Structured<double>& R
    )
{
    const int DOF = 2;
    double Re = Physics::GetInstance()->ReynoldsNumber();

    IndexRange cr = block.CellRange();
    const Structured<double>& U = block.U();
    const Structured<double>& Sxi   = block.Sxi();
    const Structured<double>& Seta  = block.Seta();
    const Structured<double>& Szeta = block.Szeta();
    const Structured<double>& Vol = block.Vol();
    const Structured<double>& MuK = block.MuK();
    const Structured<double>& TurMuK = block.TurMuK();

    ConvectiveFluxEvaluator fluxEval;
    DiffusiveFluxEvaluator diffusionEval(Re);
    double H[2], Hv[2];

    R = 0.0;

    // Convective & diffusive fluxes

    // Xi flux
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I - 1; i <= cr.End.I; ++i)
            {
                double *UL, *UR, *UTL, *UTR;
                double *RL, *RR;
                UL = U(i, j, k);
                UR = U(i + 1, j, k);
                UTL = UT(i, j, k);
                UTR = UT(i + 1, j, k);
                Vector3 Sn = Vector3(Sxi(i, j, k));
                fluxEval.EvaluateFlux(H, UL, UR, UTL, UTR, Sn);

                AuxiliaryCell cell(I, IndexIJK(i, j, k), Sxi, Seta, Szeta, Vol, cr);
                diffusionEval.EvaluateFlux(Hv, cell, UT, U, MuK, TurMuK, Sn);

                if (i > cr.Start.I - 1)
                {
                    RL = R(i, j, k);
                    for (int l = 0; l < DOF; ++l)
                    {
                        RL[l] += (H[l] - Hv[l]);
                    }
                }
                if (i < cr.End.I)
                {
                    RR = R(i + 1, j, k);
                    for (int l = 0; l < DOF; ++l)
                    {
                        RR[l] -= (H[l] - Hv[l]);
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
                double *UL, *UR, *UTL, *UTR;
                double *RL, *RR;
                UL = U(i, j, k);
                UR = U(i, j + 1, k);
                UTL = UT(i, j, k);
                UTR = UT(i, j + 1, k);
                Vector3 Sn = Vector3(Seta(i, j, k));
                fluxEval.EvaluateFlux(H, UL, UR, UTL, UTR, Sn);

                AuxiliaryCell cell(J, IndexIJK(i, j, k), Sxi, Seta, Szeta, Vol, cr);
                diffusionEval.EvaluateFlux(Hv, cell, UT, U, MuK, TurMuK, Sn);

                if (j > cr.Start.J - 1)
                {
                    RL = R(i, j, k);
                    for (int l = 0; l < DOF; ++l)
                    {
                        RL[l] += (H[l] - Hv[l]);
                    }
                }
                if (j < cr.End.J)
                {
                    RR = R(i, j + 1, k);
                    for (int l = 0; l < DOF; ++l)
                    {
                        RR[l] -= (H[l] - Hv[l]);
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
                double *UL, *UR, *UTL, *UTR;
                double *RL, *RR;
                UL = U(i, j, k);
                UR = U(i, j, k + 1);
                UTL = UT(i, j, k);
                UTR = UT(i, j, k + 1);
                Vector3 Sn = Vector3(Szeta(i, j, k));
                fluxEval.EvaluateFlux(H, UL, UR, UTL, UTR, Sn);

                AuxiliaryCell cell(K, IndexIJK(i, j, k), Sxi, Seta, Szeta, Vol, cr);
                diffusionEval.EvaluateFlux(Hv, cell, UT, U, MuK, TurMuK, Sn);

                if (k > cr.Start.K - 1)
                {
                    RL = R(i, j, k);
                    for (int l = 0; l < DOF; ++l)
                    {
                        RL[l] += (H[l] - Hv[l]);
                    }
                }
                if (k < cr.End.K)
                {
                    RR = R(i, j, k + 1);
                    for (int l = 0; l < DOF; ++l)
                    {
                        RR[l] -= (H[l] - Hv[l]);
                    }
                }
            }
        }
    }

    // Production & Dissipation
    VelocityAdaptor QAdaptor(U);
    PrimaryCellAdaptor<3, VelocityAdaptor> CellAdaptor(QAdaptor, block);
    GradientEvaluator<3, PrimaryCellAdaptor<3, VelocityAdaptor> > gradEval(CellAdaptor);

    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double mut;
                double dUdX[3][3];
                gradEval.Evaluate(dUdX, IndexIJK(i, j, k));
                mut = TurMuK(i, j, k)[0];

                double rho, tke, omega;
                rho = U(i, j, k)[0];
                tke = UT(i, j, k)[0] / rho;
                omega = UT(i, j, k)[1] / rho;

                double P[2];
                P[0] = P[1] = 0.0;
                for (int ii = 0; ii < 3; ++ii)
                {
                    for (int jj = 0; jj < 3; ++jj)
                    {
                        double Sij = 0.5 * (dUdX[ii][jj] + dUdX[jj][ii]);
                        double two_S_dUdX = 2.0 * Sij * dUdX[ii][jj];
                        double tau_dUdX = mut * two_S_dUdX;

                        P[0] += tau_dUdX;
                        P[1] += ALPHA_OMEGA * rho * two_S_dUdX;
                    }
                }
                P[0] /= Re;
                P[1] /= Re;

                double D[2];
                D[0] = BETA_K * rho * tke * omega;
                D[1] = BETA_OMEGA * rho * omega * omega;
                D[0] *= Re;
                D[1] *= Re;

                // Limit the production term
                if (P[0] < 0.0 || D[0] < 0.0 || P[1] < 0.0 || D[1] < 0.0)
                {
                    std::cout << "Negative P or D at " << IndexIJK(i, j, k) << ": " << P[0] << ", " << P[1] << ", " << D[0] << ", " << D[1] << std::endl;
                }
                P[0] = std::max(0.0, P[0]);
                P[1] = std::max(0.0, P[1]);
                D[0] = std::max(0.0, D[0]);
                D[1] = std::max(0.0, D[1]);

                P[0] = std::min(P[0], 20.0 * D[0]);

                double* r = R(i, j, k);
                double vol = *Vol(i, j, k);
                r[0] -= vol * (P[0] - D[0]);
                r[1] -= vol * (P[1] - D[1]);
            }
        }
    }

    int dof = UT.DOF();
    double square[dof], maxSquare[dof + 1];
    std::vector<IndexIJK> maxindices;
    R.ReduceSquared(square, maxSquare, maxindices, cr);
    Residual::BlockIndices maxlocs(dof + 1);
    for (int i = 0; i < dof + 1; ++i)
    {
        maxlocs[i] = Residual::BlockIndex(block.ID(), maxindices[i]);
    }
    return Residual(dof, square, maxSquare, maxlocs, cr.Count());
}

TurbulenceModel::TurbulenceModel()
{
}

void
TurbulenceModel::ApplyBCs(Block& block) const
{
    block.ApplyTurbBCs();
}

#include <cstdlib>

double
TurbulenceModel::ComputeTimeStep(
    const Block& block,
    const Structured<double>& UT,
    Structured<double>& DT,
    double cfl,
    bool localTimeStepping
    )
{
    double Re = Physics::GetInstance()->ReynoldsNumber();
    const Structured<double>& Sxi = block.Sxi();
    const Structured<double>& Seta = block.Seta();
    const Structured<double>& Szeta = block.Szeta();
    const Structured<double>& Vol = block.Vol();
    const Structured<double>& U = block.U();
    const Structured<double>& MuK = block.MuK();
    const Structured<double>& TurMuK = block.TurMuK();

    IndexRange cr = block.CellRange();

    double dtmin = 1.0e20;
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                Vector3 sxi, seta, szeta;
                sxi = 0.5 * (Vector3(Sxi(i - 1, j, k)) + Vector3(Sxi(i, j, k)));
                seta = 0.5 * (Vector3(Seta(i, j - 1, k)) + Vector3(Seta(i, j, k)));
                szeta = 0.5 * (Vector3(Szeta(i, j, k - 1)) + Vector3(Szeta(i, j, k)));

                double* uu = U(i, j, k);
                Vector3 V(uu[1] / uu[0], uu[2] / uu[0], uu[3] / uu[0]);
                double Vxi, Veta, Vzeta;
                Vxi = std::abs(dot_product(V, sxi));
                Veta = std::abs(dot_product(V, seta));
                Vzeta = std::abs(dot_product(V, szeta));

                double rho, mu, turmu, lv, vol, sxiAbs, setaAbs, szetaAbs, dxxi, dxeta, dxzeta;
                rho = uu[0];
                mu = MuK(i, j, k)[0];
                turmu = TurMuK(i, j, k)[0];
                lv = (mu + 2.0 * turmu) / rho / Re;
                vol = *Vol(i, j, k);
                sxiAbs   = sxi.Mag();   dxxi   = vol / sxiAbs;
                setaAbs  = seta.Mag();  dxeta  = vol / setaAbs;
                szetaAbs = szeta.Mag(); dxzeta = vol / szetaAbs;

                double lambda;
                //lambda = std::max(Vxi, std::max(Veta, Vzeta));
                lambda = Vxi + Veta + Vzeta + sxiAbs * lv / dxxi + setaAbs * lv / dxeta + szetaAbs * lv / dxzeta;
                double dt;
                dt = cfl * (*Vol(i, j, k)) / lambda;
                *DT(i, j, k) = dt;
                dtmin = std::min(dt, dtmin);
            }
        }
    }

    if (!localTimeStepping)
    {
        DT = dtmin;
    }

    return dtmin;
}

double
TurbulenceModel::ScalarCoeff(
    const Block& block, const Structured<double>& UT,
    const IndexIJK& Ii, const IndexIJK& Ij,
    double* Sij, double SijSign, double SijAbs, const Structured<double>& Radius
    ) const
{
    const Structured<double>& Vol = block.Vol();
    const Structured<double>& U = block.U();
    const Structured<double>& MuK = block.MuK();
    const Structured<double>& TurMuK = block.TurMuK();
    double *Ui, *Uj;
    double lambda;

    Ui = U(Ii);
    Uj = U(Ij);

#if 1
    double Vi, Vj, Vij;
    Vi = SijSign * (Sij[0] * Ui[1] + Sij[1] * Ui[2] + Sij[2] * Ui[3]) / Ui[0];
    Vj = SijSign * (Sij[0] * Uj[1] + Sij[1] * Uj[2] + Sij[2] * Uj[3]) / Uj[0];
    Vij = 0.5 * (Vi + Vj);

    lambda = std::abs(Vij);
#else
    Vector3 Sn(Sij);
    Sn = SijSign * Sn;

    Vector3 Vi(Ui[1] / Ui[0], Ui[2] / Ui[0], Ui[3] / Ui[0]);
    Vector3 Vj(Uj[1] / Uj[0], Uj[2] / Uj[0], Uj[3] / Uj[0]);
    Vector3 Vij = 0.5 * (Vi + Vj);
    double Vnij = dot_product(Vij, Sij);

    lambda = std::abs(Vnij);
#endif

    double rho, mu, mut, dx;
    rho = 0.5 * (Ui[0] + Uj[0]);
    mu = 0.5 * (MuK(Ii)[0] + MuK(Ij)[0]);
    mut = 0.5 * (TurMuK(Ii)[0] + TurMuK(Ij)[0]);
    dx = *Vol(Ii) / SijAbs;

    double Re = Physics::GetInstance()->ReynoldsNumber();
    double lambdaV = SijAbs * (mu + 2.0 * mut) / rho / dx / Re;

    return lambda + lambdaV;
}

inline
void Flux(double* H, double* U, double* UT, double* Sij, double SijSign)
{
    double Vij;

    Vij = SijSign * (Sij[0] * U[1] + Sij[1] * U[2] + Sij[2] * U[3]) / U[0];
    H[0] = Vij * UT[0];
    H[1] = Vij * UT[1];
}

void
TurbulenceModel::JacobianDU(
    double* JacDUT, const Block& block, const IndexIJK& i,
    double* UT, double* dUT, double* Sij, double SijSign, double SijAbs, double* Radius
    ) const
{
    double eps = 1.0e-8;
    double H[2], H2[2], UT2[2];
    double* U = block.U()(i);

    for (int l = 0; l < 2; ++l)
        UT2[l] = UT[l] + eps * dUT[l];

    Flux(H, U, UT, Sij, SijSign);
    Flux(H2, U, UT2, Sij, SijSign);

    for (int l = 0; l < 2; ++l)
        JacDUT[l] = (H2[l] - H[l]) / eps;
}

void
TurbulenceModel::SetUT(Structured<double>& UT, const Structured<double>& U, const IndexRange& range, const TurbulenceSpec& turbSpec) const
{
    double tke = turbSpec.TKE_Nondimensional();
    double omega = turbSpec.Omega_Nondimensional();
    for (IndexIterator itor(range); !itor.IsEnd(); itor.Advance())
    {
        IndexIJK i = itor.Index();
        const double* u = U(i);
        double* ut = UT(i);
        ut[0] = u[0] * tke;
        ut[1] = u[0] * omega;
    }
}

