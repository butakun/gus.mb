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

inline
RoeFluxEvaluator::RoeFluxEvaluator(double epsEntropyFix)
:   mEps(epsEntropyFix)
{
}

#ifdef BLAZEK_ROE
inline
void
RoeFluxEvaluator::EvaluateFlux(
    double* UL,
    double* UR,
    double* Sn,
    double Gamma,
    double* RadiusL,
    double* RadiusR,
    double* H
    ) const
{
    double kx, ky, kz, kabs, dkabs, nx, ny, nz;
    kx    = Sn[0];
    ky    = Sn[1];
    kz    = Sn[2];
    kabs  = std::sqrt(kx * kx + ky * ky + kz * kz);
    dkabs = 1.0 / kabs;
    nx    = kx * dkabs;
    ny    = ky * dkabs;
    nz    = kz * dkabs;

    // Left state
    double rhoL, uL, vL, wL, keL, romegaSqL, kerL, rhoeL, gammaL, pL, hL, htL, VnL, rhoVL;
    rhoL = UL[0];
    uL = UL[1] / rhoL;
    vL = UL[2] / rhoL;
    wL = UL[3] / rhoL;
    keL = 0.5 * (uL * uL + vL * vL + wL * wL);
    romegaSqL = RadiusL[3]; // (Radius * omega)**2
    kerL = 0.5 * romegaSqL;
    rhoeL = UL[4] - rhoL * keL + rhoL * kerL;
    gammaL = Gamma;
    pL = (gammaL - 1.0) * rhoeL;
    hL = rhoeL * gammaL / rhoL;
    htL = hL + keL - kerL; // rothalpy
    VnL = nx * uL + ny * vL + nz * wL;
    rhoVL = rhoL * (kx * uL + ky * vL + kz * wL);

    // Right state
    double rhoR, uR, vR, wR, keR, romegaSqR, kerR, rhoeR, gammaR, pR, hR, htR, VnR, rhoVR;
    rhoR = UR[0];
    uR = UR[1] / rhoR;
    vR = UR[2] / rhoR;
    wR = UR[3] / rhoR;
    keR = 0.5 * (uR * uR + vR * vR + wR * wR);
    romegaSqR = RadiusR[3]; // (Radius * omega)**2
    kerR = 0.5 * romegaSqR;
    rhoeR = UR[4] - rhoR * keR + rhoR * kerR;
    gammaR = Gamma;
    pR = (gammaR - 1.0) * rhoeR;
    hR = rhoeR * gammaR / rhoR;
    htR = hR + keR - kerR; // rothalpy
    VnR = nx * uR + ny * vR + nz * wR;
    rhoVR = rhoR * (kx * uR + ky * vR + kz * wR);

    // Roe averated state
    double R, DRP1, rho, u, v, w, ht, gamma, Vn, qsq, ker;
    R     = std::sqrt(rhoR / rhoL);
    DRP1  = 1.0 / (R + 1.0);
    rho   = R * rhoL;
    u     = (R * uR + uL) * DRP1;
    v     = (R * vR + vL) * DRP1;
    w     = (R * wR + wL) * DRP1;
    ht    = (R * htR + htL) * DRP1;
    gamma = 0.5 * (gammaL + gammaR);
    Vn    = nx * u + ny * v + nz * w;
    qsq   = u * u + v * v + w * w;
    ker   = 0.5 * (romegaSqL + romegaSqR);

    // Primitive variable differences
    double drho, du, dv, dw, dp, dVn;
    drho = rhoR - rhoL;
    du   = uR - uL;
    dv   = vR - vL;
    dw   = wR - wL;
    dp   = pR - pL;
    dVn  = VnR - VnL;

    // Misc related variables
    double gm1, dgm1, ke, csq, c, dc, dcsq, rhoc, drhoc, pLR;
    gm1   = gamma - 1.0;
    dgm1  = 1.0 / gm1;
    ke    = 0.5 * qsq;
    csq   = gm1 * (ht - ke + ker);
    c     = std::sqrt(csq);
    dc    = 1.0 / c;
    dcsq  = 1.0 / csq;
    rhoc  = rho * c;
    drhoc = 1.0 / rhoc;
    pLR   = pL + pR;

    // Eigenvalues
    double VnAbs, VpCabs, VmCabs;
    VnAbs  = std::max(std::abs(Vn), mEps);
    VpCabs = std::max(std::abs(Vn + c), mEps);
    VmCabs = std::max(std::abs(Vn - c), mEps);

    // Roe correction term
    double dF[5];
    double C1, C234, C5;
    C1   = kabs * VmCabs * (dp - rhoc * dVn) * 0.5 * dcsq;
    C234 = kabs * VnAbs * (drho - dp * dcsq);
    C5   = kabs * VpCabs * (dp + rhoc * dVn) * 0.5 * dcsq;

    dF[0] = C1;
    dF[1] = C1 * (u - nx * c);
    dF[2] = C1 * (u - ny * c);
    dF[3] = C1 * (u - nz * c);
    dF[4] = C1 * (ht - c * Vn);

    dF[0] += C234;
    dF[1] += C234 * u + kabs * VnAbs * rho * (du - dVn * nx);
    dF[2] += C234 * v + kabs * VnAbs * rho * (dv - dVn * ny);
    dF[3] += C234 * w + kabs * VnAbs * rho * (dw - dVn * nz);
    dF[4] += C234 * ke + kabs * VnAbs * rho * (u * du + v * dv + w * dw - Vn * dVn);

    dF[0] += C5;
    dF[1] += C5 * (u + c * nx);
    dF[2] += C5 * (v + c * ny);
    dF[3] += C5 * (w + c * nz);
    dF[4] += C5 * (ht + c * Vn);

    // Roe flux
    H[0] = 0.5 * (rhoVL + rhoVR                     ) - 0.5 * dF[0];
    H[1] = 0.5 * (rhoVL * uL + rhoVR * uR + kx * pLR) - 0.5 * dF[1];
    H[2] = 0.5 * (rhoVL * vL + rhoVR * vR + ky * pLR) - 0.5 * dF[2];
    H[3] = 0.5 * (rhoVL * wL + rhoVR * wR + kz * pLR) - 0.5 * dF[3];
    H[4] = 0.5 * (htL * rhoVL + htR * rhoVR         ) - 0.5 * dF[4];
}

#else

inline
void
RoeFluxEvaluator::EvaluateFlux(
    double* UL,
    double* UR,
    double* Sn,
    double Gamma,
    double* RadiusL,
    double* RadiusR,
    double* H
    ) const
{
    double kx, ky, kz, kabs, dkabs, nx, ny, nz;
    kx    = Sn[0];
    ky    = Sn[1];
    kz    = Sn[2];
    kabs  = std::sqrt(kx * kx + ky * ky + kz * kz);
    dkabs = 1.0 / kabs;
    nx    = kx * dkabs;
    ny    = ky * dkabs;
    nz    = kz * dkabs;

    // Left state
    double rhoL, uL, vL, wL, keL, romegaSqL, kerL, rhoeL, gammaL, pL, hL, htL, rhoVL;
    rhoL = UL[0];
    uL = UL[1] / rhoL;
    vL = UL[2] / rhoL;
    wL = UL[3] / rhoL;
    keL = 0.5 * (uL * uL + vL * vL + wL * wL);
    romegaSqL = RadiusL[3];
    kerL = 0.5 * romegaSqL;
    rhoeL = UL[4] - rhoL * keL + rhoL * kerL;
    gammaL = Gamma;
    pL = (gammaL - 1.0) * rhoeL;
    hL = gammaL * rhoeL / rhoL;
    htL = hL + keL - kerL; // rothalpy
    rhoVL = rhoL * (kx * uL + ky * vL + kz * wL);

    // Right state
    double rhoR, uR, vR, wR, keR, romegaSqR, kerR, rhoeR, gammaR, pR, hR, htR, rhoVR;
    rhoR = UR[0];
    uR = UR[1] / rhoR;
    vR = UR[2] / rhoR;
    wR = UR[3] / rhoR;
    keR = 0.5 * (uR * uR + vR * vR + wR * wR);
    romegaSqR = RadiusR[3];
    kerR = 0.5 * romegaSqR;
    rhoeR = UR[4] - rhoR * keR + rhoR * kerR;
    gammaR = Gamma;
    pR = (gammaR - 1.0) * rhoeR;
    hR = gammaR * rhoeR / rhoR;
    htR = hR + keR - kerR; // rothalpy
    rhoVR = rhoR * (kx * uR + ky * vR + kz * wR);

    // Roe averated state
    double R, DRP1, rho, u, v, w, ht, gamma;
    R     = std::sqrt(rhoR / rhoL);
    DRP1  = 1.0 / (R + 1.0);
    rho   = R * rhoL;
    u     = (R * uR + uL) * DRP1;
    v     = (R * vR + vL) * DRP1;
    w     = (R * wR + wL) * DRP1;
    ht    = (R * htR + htL) * DRP1;
    gamma = 0.5 * (gammaL + gammaR);

    // Primitive variable differences
    double drho, du, dv, dw, dp;
    drho = rhoR - rhoL;
    du   = uR - uL;
    dv   = vR - vL;
    dw   = wR - wL;
    dp   = pR - pL;

    // Misc related variables
    double gm1, dgm1, ke, ker, csq, c, dc, dcsq, drhoc, Vk, VkHat, dVkHat, ck;
    gm1   = gamma - 1.0;
    dgm1  = 1.0 / gm1;
    ke    = 0.5 * (u * u + v * v + w * w);
    ker   = 0.5 * (kerL + kerR);
    csq   = gm1 * (ht - ke + ker);
    c     = std::sqrt(csq);
    dc    = 1.0 / c;
    dcsq  = 1.0 / csq;
    drhoc = 1.0 / (rho * c);

    Vk     = kx * u + ky * v + kz * w;
    VkHat  = Vk * dkabs;
    dVkHat = nx * du + ny * dv + nz * dw;
    ck     = kabs * c;

    // Harten's entropy fix parameter scaled to the metric
    double epsk;
    epsk   = mEps * kabs;

    // Eigenvalues
    double VkAbs, VpCabs, VmCabs;
    VkAbs  = std::max(std::abs(Vk), epsk);
    VpCabs = std::max(std::abs(Vk + ck), epsk);
    VmCabs = std::max(std::abs(Vk - ck), epsk);

    // Another misc related variables
    double pLR, dpdrhoc, rhod2c, dpdcsq, rhoVkAbs, DDp, DDm, B1, B2, B3, B4, B5, B6, B7;
    pLR      = pL + pR;
    dpdrhoc  = dp * drhoc;
    rhod2c   = 0.5 * rho * dc;
    dpdcsq   = dp * dcsq;
    rhoVkAbs = rho * VkAbs;
    DDp      = VpCabs * (dpdrhoc + dVkHat);
    DDm      = VmCabs * (dpdrhoc - dVkHat);
    B1       = VkAbs * (drho - dpdcsq);
    B2       = rhod2c * (DDp + DDm);
    B3       = B1 + B2;
    B4       = 0.5 * rho * (DDp - DDm);
    B5       = rhoVkAbs * (du - nx * dVkHat);
    B6       = rhoVkAbs * (dv - ny * dVkHat);
    B7       = rhoVkAbs * (dw - nz * dVkHat);

    // Roe correction term
    double dF1, dF2, dF3, dF4, dF5;
    dF1 = B3;
    dF2 = u * B3 + B5 + nx * B4;
    dF3 = v * B3 + B6 + ny * B4;
    dF4 = w * B3 + B7 + nz * B4;
    dF5 = ht * B3 - csq * dgm1 * B1 + VkHat * B4 + u * B5 + v * B6 + w * B7;

    // Roe flux
    H[0] = 0.5 * (rhoVL + rhoVR                      - dF1);
    H[1] = 0.5 * (rhoVL * uL + rhoVR * uR + kx * pLR - dF2);
    H[2] = 0.5 * (rhoVL * vL + rhoVR * vR + ky * pLR - dF3);
    H[3] = 0.5 * (rhoVL * wL + rhoVR * wR + kz * pLR - dF4);
    H[4] = 0.5 * (htL * rhoVL + htR * rhoVR          - dF5);
}
#endif

