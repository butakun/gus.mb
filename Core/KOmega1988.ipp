// $Id: KOmega1988.ipp 284 2013-06-14 03:22:13Z kato $

#include "AuxiliaryCell.h"
#include "GradientEvaluator.h"
#include "Physics.h"

using namespace KOmega1988;

inline
void
KOmega1988::ConvectiveFluxEvaluator::EvaluateFlux(
    double* H, double* UL, double* UR, double* UTL, double* UTR, const Vector3& Sn
    )
{
    Vector3 VelL(UL[1] / UL[0], UL[2] / UL[0], UL[3] / UL[0]);
    Vector3 VelR(UR[1] / UR[0], UR[2] / UR[0], UR[3] / UR[0]);
    //Vector3 VelLR = 0.5 * (VelL + VelR);
    //Vector3 n = Sn.Normalized();
    double SnAbs = Sn.Mag();

    double VnL = dot_product(VelL, Sn);
    double VnR = dot_product(VelR, Sn);
#if 1
    double VnLR = 0.5 * (VnL + VnR);
#else
    double VnLR = dot_product(VelLR, Sn);
#endif

    double VnLRAbs = std::abs(VnLR);
    //double VnSign = VnLR >= 0.0 ? 1.0 : -1.0;

    double HL[2], HR[2];

    HL[0] = VnL * UTL[0];
    HL[1] = VnL * UTL[1];
    HR[0] = VnR * UTR[0];
    HR[1] = VnR * UTR[1];

    double eps = 0.0 * 0.02 * SnAbs;

    H[0] = 0.5 * (HR[0] + HL[0]) - 0.5 * (VnLRAbs + eps) * (UTR[0] - UTL[0]);
    H[1] = 0.5 * (HR[1] + HL[1]) - 0.5 * (VnLRAbs + eps) * (UTR[1] - UTL[1]);
}

inline
void
KOmega1988::DiffusiveFluxEvaluator::EvaluateFlux(
    double* H,
    const AuxiliaryCell& cell,
    const Structured<double>& UT,
    const Structured<double>& U,
    const Structured<double>& MuK,
    const Structured<double>& TurMuK,
    const Vector3& Sn
    )
{
    double mu, mut;
    IndexIJK i1, i2;

    i1 = cell.IndexU1M();
    i2 = cell.IndexU1P();
    mu = 0.5 * (MuK(i1)[0] + MuK(i2)[0]);
    mut = 0.5 * (TurMuK(i1)[0] + TurMuK(i2)[0]);

    KOmegaAdaptor QAdaptor(U, UT);
    AuxiliaryCellAdaptor<2, KOmegaAdaptor> CellAdaptor(cell, QAdaptor);
    GradientEvaluator<2, AuxiliaryCellAdaptor<2, KOmegaAdaptor> > gradEval(CellAdaptor);

    double dUTdX[2][3];
    gradEval.Evaluate(dUTdX, i1);

    double dkdx, dkdy, dkdz, dodx, dody, dodz;
    dkdx = dUTdX[0][0]; dkdy = dUTdX[0][1]; dkdz = dUTdX[0][2];
    dodx = dUTdX[1][0]; dody = dUTdX[1][1]; dodz = dUTdX[1][2];

    H[0] = 1.0 / mRe * (mu +     SIGMA_K * mut) * (Sn.X() * dkdx + Sn.Y() * dkdy + Sn.Z() * dkdz);
    H[1] = 1.0 / mRe * (mu + SIGMA_OMEGA * mut) * (Sn.X() * dodx + Sn.Y() * dody + Sn.Z() * dodz);
}

#if SCALAR_DIAG
inline
double
TurbulenceModel::DiagonalAddition(const Block& block, const Structured<double>& UT, const IndexIJK& Ii) const
{
    double Re = Physics::GetInstance()->ReynoldsNumber();

    const Structured<double>& U = block.U();
    const Structured<double>& Vol = block.Vol();

    return (0.9 + 0.0) * UT(Ii)[1] / U(Ii)[0] * Re * (*Vol(Ii));
}
#else
inline
void
TurbulenceModel::DiagonalAddition(double* D, const Block& block, const Structured<double>& UT, const IndexIJK& Ii) const
{
    double Re = Physics::GetInstance()->ReynoldsNumber(); // FIXME: we could cache this in the ctor?

    const Structured<double>& U = block.U();
    const Structured<double>& Vol = block.Vol();

    double vol = *Vol(Ii);
    double omega = UT(Ii)[1] / U(Ii)[0];

    D[0] += vol * Re * BETA_K * omega;
    D[1] += vol * Re * 2.0 * BETA_OMEGA * omega;
}
#endif

#include "Physics.h"

inline
double
KOmega1988::TurbulenceModel::MuT(const double* U, const double* UT) const
{
    //double Re = Physics::GetInstance()->ReynoldsNumber();
    double rho = U[0], rhok = UT[0], rhoomega = UT[1];

    //double mut = Re * rho * rhok / rhoomega;
    rhok = std::max(0.0, rhok);
    rhoomega = std::max(1.0e-20, rhoomega);
    double mut = rho * rhok / rhoomega;
    mut = std::max(0.0, mut);

    return mut;
}

