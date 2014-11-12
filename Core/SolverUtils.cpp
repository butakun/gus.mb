// $Id: SolverUtils.cpp 129 2011-09-08 16:42:22Z kato $

#include "SolverUtils.h"
#include "Physics.h"
#include "Vector3.h"

void
SolverUtils::ComputeSpectralRadius(
    const Block& block,
    const Structured<double>& U,
    Structured<double>& Lambda
    )
{
    IndexRange cellRange = block.CellRange();

    const Structured<double>& Sxi = block.Sxi();
    const Structured<double>& Seta = block.Seta();
    const Structured<double>& Szeta = block.Szeta();
    const Structured<double>& Vol = block.Vol();
    const Structured<double>& MuK = block.MuK();
    const Structured<double>& TurMuK = block.TurMuK();

    double gamma = Physics::GetInstance()->Gamma();
    double Re = Physics::GetInstance()->ReynoldsNumber();

    for (int k = cellRange.Start.K; k <= cellRange.End.K; ++k)
    {
        for (int j = cellRange.Start.J; j <= cellRange.End.J; ++j)
        {
            for (int i = cellRange.Start.I; i <= cellRange.End.I; ++i)
            {
                Vector3 sxi   = 0.5 * (Vector3(  Sxi(i - 1, j, k)) + Vector3(  Sxi(i, j, k)));
                Vector3 seta  = 0.5 * (Vector3( Seta(i, j - 1, k)) + Vector3( Seta(i, j, k)));
                Vector3 szeta = 0.5 * (Vector3(Szeta(i, j, k - 1)) + Vector3(Szeta(i, j, k)));

                double sxiAbs   = sxi.Mag();
                double setaAbs  = seta.Mag();
                double szetaAbs = szeta.Mag();
                double vol = *Vol(i, j, k);

                double* UU = U(i, j, k);
                double rho, u, v, w, rhoe, p, c;
                rho = UU[0];
                u = UU[1] / rho;
                v = UU[2] / rho;
                w = UU[3] / rho;
                rhoe = UU[4] - 0.5 * rho * (u * u + v * v + w * w); // FIXME: rothalpy
                p = rhoe * (gamma - 1.0);
                c = std::sqrt(gamma * p / rho);

                double mu, mut, lv;
                mu  = MuK(i, j, k)[0];
                mut = TurMuK(i, j, k)[0];
                lv = (mu + mut) / rho / Re;

                double lambdaXi, lambdaEta, lambdaZeta;
                lambdaXi   = std::abs(sxi[0] * u + sxi[1] * v + sxi[2] * w) + sxiAbs * c + sxiAbs * lv * sxiAbs / vol;
                lambdaEta  = std::abs(seta[0] * u + seta[1] * v + seta[2] * w) + setaAbs * c + setaAbs * lv * setaAbs / vol;
                lambdaZeta = std::abs(szeta[0] * u + szeta[1] * v + szeta[2] * w) + szetaAbs * c + szetaAbs * lv * szetaAbs / vol;

                *Lambda(i, j, k) = lambdaXi + lambdaEta + lambdaZeta;
            }
        }
    }
}

double
SolverUtils::ComputeTimeStep(
    const Block& block,
    const Structured<double>& U,
    Structured<double>& DT,
    double cfl,
    bool localStepping
    )
{
    ComputeSpectralRadius(block, U, DT);

    IndexRange cellRange = block.CellRange();

    const Structured<double>& Vol = block.Vol();

    // Initialize the minimum.
    double dtmin = cfl * *Vol(cellRange.Start) / *DT(cellRange.Start);
    //double dtmin = 1e20;

    for (int k = cellRange.Start.K; k <= cellRange.End.K; ++k)
    {
        for (int j = cellRange.Start.J; j <= cellRange.End.J; ++j)
        {
            for (int i = cellRange.Start.I; i <= cellRange.End.I; ++i)
            {
                double vol = *Vol(i, j, k);
                double lambda = *DT(i, j, k);
                double dt = cfl * vol / lambda;
                dtmin = std::min(dtmin, dt);
                *DT(i, j, k) = dt;
            }
        }
    }

    if (!localStepping)
    {
        for (int k = cellRange.Start.K; k <= cellRange.End.K; ++k)
        {
            for (int j = cellRange.Start.J; j <= cellRange.End.J; ++j)
            {
                for (int i = cellRange.Start.I; i <= cellRange.End.I; ++i)
                {
                    *DT(i, j, k) = dtmin;
                }
            }
        }
    }

    return dtmin;
}

