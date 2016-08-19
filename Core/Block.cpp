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

#include "Communicator.h"
#include "Block.h"
#include "SolverFunctions.h"
#include "Vector3.h"
#include "BC.h"
#include "Connectivity1to1.h"
#include "Physics.h"
#include "Roster.h"
#include "RigidBodyMotion.h"
#include <cassert>

Block*
Block::New(int id, const IndexRange& meshRange, bool unsteady)
{
    int temporalstore = unsteady ? 2 : 0;
    Block* block = new Block(id, meshRange, temporalstore, 2);

    Roster::GetInstance()->RegisterBlock(Communicator::GetInstance()->MyRank(), block);

    return block;
}

Block::Block(int id, const IndexRange& meshRange, int temporalstore, int nGhosts)
:   VirtualBlock(id, meshRange),
    mNGhostLayers(nGhosts),
    mDim(meshRange.Start - IndexIJK(nGhosts) + IndexIJK(1), meshRange.End + IndexIJK(nGhosts)),
    mXYZ(3, meshRange), mU(5, mDim), mUT(2, mDim),
    mSxi(3, mDim), mSeta(3, mDim), mSzeta(3, mDim), mVol(1, mDim), mRadius(4, mDim),
    mMuK(2, mDim), mTurMuK(2, mDim), mMask(1, mDim)
    //mRotating(false), mAxisOrigin(0.0, 0.0, 0.0), mAngularVelocity(1.0, 0.0, 0.0)
{
    assert(temporalstore <= 2);
    if (temporalstore > 1)
    {
        mU2.Allocate(mU.DOF(), mDim);
        mU3.Allocate(mU.DOF(), mDim);
        mUT2.Allocate(mUT.DOF(), mDim);
        mUT3.Allocate(mUT.DOF(), mDim);
    }

    assert(meshRange.Start == IndexIJK(0, 0, 0));
}

Block::~Block()
{
    delete[] mXYZ.Data;
    delete[] mU.Data;
    delete[] mU2.Data;
    delete[] mU3.Data;
    delete[] mUT.Data;
    delete[] mUT2.Data;
    delete[] mUT3.Data;
    delete[] mSxi.Data;
    delete[] mSeta.Data;
    delete[] mSzeta.Data;
    delete[] mVol.Data;
    delete[] mRadius.Data;
    delete[] mMuK.Data;
    delete[] mTurMuK.Data;
    delete[] mMask.Data;
}

void
Block::ComputeMetrics()
{
    IndexRange mr = MeshRange();
    IndexRange cr = CellRange();
    Structured<double> xyz = XYZ();
    Structured<double> sxi = Sxi();
    Structured<double> seta = Seta();
    Structured<double> szeta = Szeta();
    Structured<double> vol = Vol();
    //Structured<double> rad = Radius();

    // Sxi
    for (int i = mr.Start.I; i <= mr.End.I; ++i)
    {
        int i0 = i;
        for (int j = mr.Start.J + 1; j <= mr.End.J; ++j)
        {
            int j0 = j - 1, j1 = j;
            for (int k = mr.Start.K + 1; k <= mr.End.K; ++k)
            {
                int k0 = k - 1, k1 = k;
                Vector3 d1(xyz(i0, j0, k0), xyz(i0, j1, k1));
                Vector3 d2(xyz(i0, j1, k0), xyz(i0, j0, k1));
                Vector3 s = 0.5 * cross_product(d1, d2);
                double* sn = sxi(i, j, k);
                sn[0] = s[0];
                sn[1] = s[1];
                sn[2] = s[2];
            }
        }
    }

    // Seta
    for (int i = mr.Start.I + 1; i <= mr.End.I; ++i)
    {
        int i0 = i - 1, i1 = i;
        for (int j = mr.Start.J; j <= mr.End.J; ++j)
        {
            int j0 = j;
            for (int k = mr.Start.K + 1; k <= mr.End.K; ++k)
            {
                int k0 = k - 1, k1 = k;
                Vector3 d1(xyz(i0, j0, k0), xyz(i1, j0, k1));
                Vector3 d2(xyz(i0, j0, k1), xyz(i1, j0, k0));
                Vector3 s = 0.5 * cross_product(d1, d2);
                double* sn = seta(i, j, k);
                sn[0] = s[0];
                sn[1] = s[1];
                sn[2] = s[2];
            }
        }
    }

    // Szeta
    for (int i = mr.Start.I + 1; i <= mr.End.I; ++i)
    {
        int i0 = i - 1, i1 = i;
        for (int j = mr.Start.J + 1; j <= mr.End.J; ++j)
        {
            int j0 = j - 1, j1 = j;
            for (int k = mr.Start.K; k <= mr.End.K; ++k)
            {
                int k0 = k;
                Vector3 d1(xyz(i0, j0, k0), xyz(i1, j1, k0));
                Vector3 d2(xyz(i1, j0, k0), xyz(i0, j1, k0));
                Vector3 s = 0.5 * cross_product(d1, d2);
                double* sn = szeta(i, j, k);
                sn[0] = s[0];
                sn[1] = s[1];
                sn[2] = s[2];
            }
        }
    }

    // Volume
    for (int i = cr.Start.I; i <= cr.End.I; ++i)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int k = cr.Start.K; k <= cr.End.K; ++k)
            {
                Vector3 s_xi   = sxi(i - 1, j, k);
                Vector3 s_eta  = seta(i, j - 1, k);
                Vector3 s_zeta = szeta(i, j, k - 1);
                Vector3 ss = (s_xi + s_eta + s_zeta) / 3.0;
                Vector3 d(xyz(i - 1, j - 1, k - 1), xyz(i, j, k));
                double* v = vol(i, j, k);
                *v = dot_product(ss, d);
                if (*v <= 0.0)
                {
                    std::cout << "Negative or zero volume at " << IndexIJK(i, j, k) << std::endl;
                    assert(false);
                }
            }
        }
    }

    // Radial positions of the cells
    if (IsRotating())
    {
        double vRef = Physics::GetInstance()->VRef();
        RotationalMotion* rbm = dynamic_cast<RotationalMotion*>(GetRigidBodyMotion());
        assert(rbm != NULL);
        Vector3 omega = rbm->AngularVelocity();
        for (int k = cr.Start.K; k <= cr.End.K; ++k)
        {
            for (int j = cr.Start.J; j <= cr.End.J; ++j)
            {
                for (int i = cr.Start.I; i <= cr.End.I; ++i)
                {
                    Vector3 p1(xyz(i - 1, j - 1, k - 1));
                    Vector3 p2(xyz(i,     j - 1, k - 1));
                    Vector3 p3(xyz(i,     j    , k - 1));
                    Vector3 p4(xyz(i - 1, j    , k - 1));
                    Vector3 p5(xyz(i - 1, j - 1, k    ));
                    Vector3 p6(xyz(i,     j - 1, k    ));
                    Vector3 p7(xyz(i,     j    , k    ));
                    Vector3 p8(xyz(i - 1, j    , k    ));
                    Vector3 pC = 0.125 * (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8);
                    Vector3 vR = rbm->RadialVector(pC);
                    Vector3 omegaR = cross_product(omega, vR) / vRef;
                    mRadius(i, j, k)[0] = pC.X();
                    mRadius(i, j, k)[1] = pC.Y();
                    mRadius(i, j, k)[2] = pC.Z();
                    mRadius(i, j, k)[3] = omegaR.MagSq();
                }
            }
        }
    }
    else
    {
        mRadius = 0.0;
    }

    // Copy some metrics to the rind cells. FIXME: we should read the mesh and compute the real values.
    for (int j = cr.Start.J; j <= cr.End.J; ++j)
    {
        for (int k = cr.Start.K; k <= cr.End.K; ++k)
        {
            for (int g = 1; g <= GhostLayers(); ++g)
            {
                int i;
                i = cr.Start.I - g;
                vol(i, j, k)[0] = vol(i + 1, j, k)[0];
                for (int l = 0; l < 4; ++l)
                {
                    mRadius(i, j, k)[l] = mRadius(i + 1, j, k)[l];
                }
                i = cr.End.I + g;
                vol(i, j, k)[0] = vol(i - 1, j, k)[0];
                for (int l = 0; l < 4; ++l)
                {
                    mRadius(i, j, k)[l] = mRadius(i - 1, j, k)[l];
                }
            }
        }
    }
    for (int i = cr.Start.I; i <= cr.End.I; ++i)
    {
        for (int k = cr.Start.K; k <= cr.End.K; ++k)
        {
            for (int g = 1; g <= GhostLayers(); ++g)
            {
                int j;
                j = cr.Start.J - g;
                vol(i, j, k)[0] = vol(i, j + 1, k)[0];
                for (int l = 0; l < 4; ++l)
                {
                    mRadius(i, j, k)[l] = mRadius(i, j + 1, k)[l];
                }
                j = cr.End.J + g;
                vol(i, j, k)[0] = vol(i, j - 1, k)[0];
                for (int l = 0; l < 4; ++l)
                {
                    mRadius(i, j, k)[l] = mRadius(i, j - 1, k)[l];
                }
            }
        }
    }
    for (int i = cr.Start.I; i <= cr.End.I; ++i)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int g = 1; g <= GhostLayers(); ++g)
            {
                int k;
                k = cr.Start.K - g;
                vol(i, j, k)[0] = vol(i, j, k + 1)[0];
                for (int l = 0; l < 4; ++l)
                {
                    mRadius(i, j, k)[l] = mRadius(i, j, k + 1)[l];
                }
                k = cr.End.K + g;
                vol(i, j, k)[0] = vol(i, j, k - 1)[0];
                for (int l = 0; l < 4; ++l)
                {
                    mRadius(i, j, k)[l] = mRadius(i, j, k - 1)[l];
                }
            }
        }
    }
}

void
Block::ComputeTransportProperties()
{
    ComputeTransportProperties(mMuK, mU, mRadius);
}

void
Block::ComputeTransportProperties(Structured<double>& MuK, Structured<double>& U, const Structured<double>& Radius) const
{
    double gamma = Physics::GetInstance()->Gamma();

    for (int k = mDim.Start.K; k <= mDim.End.K; ++k)
    {
        for (int j = mDim.Start.J; j <= mDim.End.J; ++j)
        {
            for (int i = mDim.Start.I; i <= mDim.End.I; ++i)
            {
                double* u = U(i, j, k);
                double* muk = MuK(i, j, k);
                double p, T, romegaSq, rhoker;
                romegaSq = Radius(i, j, k)[3];
                rhoker = 0.5 * u[0] * romegaSq;
                p = (gamma - 1.0) * (u[4] - 0.5 * (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0] + rhoker);
                T = p / u[0];
                muk[0] = Physics::GetInstance()->Mu(T);
                muk[1] = muk[0];
            }
        }
    }
}

#include "KOmega1988.h"

void
Block::ComputeTurbulentTransportProperties()
{
    std::ostream& LOG = Communicator::GetInstance()->Console();

    const double max_ratio = 10000.0;

    KOmega1988::TurbulenceModel model;

    for (int k = mDim.Start.K; k <= mDim.End.K; ++k)
    {
        for (int j = mDim.Start.J; j <= mDim.End.J; ++j)
        {
            for (int i = mDim.Start.I; i <= mDim.End.I; ++i)
            {
                double* u = U()(i, j, k);
                double* ut = UT()(i, j, k);
                double* turmuk = TurMuK()(i, j, k);
                double* muk = MuK()(i, j, k);
                turmuk[0] = model.MuT(u, ut);
                turmuk[1] = turmuk[0];
                double ratio = turmuk[0] / muk[0];
                if (ratio > max_ratio)
                {
                    LOG << "TurMu / Mu = " << ratio << ", TurMu = " << turmuk[0] << "(will be " << max_ratio * muk[0] << ") at (" << ID() << ") " << IndexIJK(i, j, k) << std::endl;
                    turmuk[0] = max_ratio * muk[0];
                    turmuk[1] = turmuk[0];
                }
            }
        }
    }
}

bool
Block::CheckNegatives(const Structured<double>& U, std::vector<IndexIJK>& indices) const
{
    std::ostream& log = Communicator::GetInstance()->Console();

    IndexRange range = CellRange();

    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                double* u = U(i, j, k);
                double romegaSq = mRadius(i, j, k)[3];
                double rhoe;
                rhoe = u[4] - 0.5 * (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0] + 0.5 * u[0] * romegaSq;
                if (rhoe < 0.0 || rhoe != rhoe)
                {
                    indices.push_back(IndexIJK(i, j, k));
                }
            }
        }
    }

    bool negative = indices.size() > 0;
    if (negative)
    {
        log << "Negative quantities found in the following cells:" << std::endl;
        for (size_t i = 0; i < indices.size(); ++i)
        {
            IndexIJK ijk = indices[i];
            double* u = U(ijk);
            double* rads = mRadius(ijk);
            log << "Block " << ID() << ", " << ijk << ", " << Vector3(mXYZ(ijk)) << ", " << u[0] << ' ' << u[1] << ' ' << u[2] << ' ' << u[3] << ' ' << u[4] << ", romega2 = " << rads[3] << std::endl;
        }
        std::cout << "There are negative cells in Block " << ID() << std::endl;
    }

    return negative;
}

void
Block::RegisterBC(BC* bc)
{
    mBCs.push_back(bc);
}

void
Block::ApplyBCs()
{
    ApplyBCs(U());
}

void
Block::ApplyBCs(Structured<double> U)
{
    for (BCs::iterator i = mBCs.begin(); i != mBCs.end(); ++i)
    {
        BC* bc = *i;
        bc->Apply(*this, U);
    }

    FillCornerGhosts();
}

void
Block::ApplyTurbBCs()
{
    for (BCs::iterator i = mBCs.begin(); i != mBCs.end(); ++i)
    {
        BC* bc = *i;
        bc->ApplyTurb(*this, mUT, mU);
    }

    FillCornerGhosts();
}

#include "BCViscousWall.h"
void
Block::ApplyViscousWallBC()
{
    for (BCs::iterator i = mBCs.begin(); i != mBCs.end(); ++i)
    {
        BC* bc = *i;
        BCViscousWall* bcvw = dynamic_cast<BCViscousWall*>(bc);
        if (bcvw == NULL)
            continue;
        bcvw->ApplyTurMuK(*this);
    }
}

void
Block::FinalizeConnectivities(Structured<double>& U)
{
    // FIXME:
    for (BCs::iterator i = mBCs.begin(); i != mBCs.end(); ++i)
    {
        BC* bc = *i;
        Connectivity1to1* conn1to1 = dynamic_cast<Connectivity1to1*>(bc);
        if (conn1to1 == NULL)
        {
            continue;
        }
        conn1to1->ApplyPeriodicity(U);
    }
}

void
Block::FillCornerGhosts()
{
    FillCornerGhosts(mU);
    FillCornerGhosts(mUT);
    FillCornerGhosts(mMuK);
    FillCornerGhosts(mTurMuK);
}

void
Block::FillCornerGhosts(Structured<double>& U) const
{
    int dof = U.DOF();
    IndexRange cr = CellRange();

    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        int i, j;
        double *u, *u1, *u2;

        i = cr.Start.I - 1;
        j = cr.Start.J - 1;
        u = U(i, j, k);
        u1 = U(i + 1, j, k);
        u2 = U(i, j + 1, k);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }

        i = cr.Start.I - 1;
        j = cr.End.J + 1;
        u = U(i, j, k);
        u1 = U(i + 1, j, k);
        u2 = U(i, j - 1, k);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }

        i = cr.End.I + 1;
        j = cr.Start.J - 1;
        u = U(i, j, k);
        u1 = U(i - 1, j, k);
        u2 = U(i, j + 1, k);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }

        i = cr.End.I + 1;
        j = cr.End.J + 1;
        u = U(i, j, k);
        u1 = U(i - 1, j, k);
        u2 = U(i, j - 1, k);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }
    }

    for (int j = cr.Start.J; j <= cr.End.J; ++j)
    {
        int i, k;
        double *u, *u1, *u2;

        i = cr.Start.I - 1;
        k = cr.Start.K - 1;
        u = U(i, j, k);
        u1 = U(i + 1, j, k);
        u2 = U(i, j, k + 1);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }

        i = cr.Start.I - 1;
        k = cr.End.K + 1;
        u = U(i, j, k);
        u1 = U(i + 1, j, k);
        u2 = U(i, j, k - 1);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }

        i = cr.End.I + 1;
        k = cr.Start.K - 1;
        u = U(i, j, k);
        u1 = U(i - 1, j, k);
        u2 = U(i, j, k + 1);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }

        i = cr.End.I + 1;
        k = cr.End.K + 1;
        u = U(i, j, k);
        u1 = U(i - 1, j, k);
        u2 = U(i, j, k - 1);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }
    }

    for (int i = cr.Start.I; i <= cr.End.I; ++i)
    {
        int j, k;
        double *u, *u1, *u2;

        j = cr.Start.J - 1;
        k = cr.Start.K - 1;
        u = U(i, j, k);
        u1 = U(i, j + 1, k);
        u2 = U(i, j, k + 1);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }

        j = cr.Start.J - 1;
        k = cr.End.K + 1;
        u = U(i, j, k);
        u1 = U(i, j + 1, k);
        u2 = U(i, j, k - 1);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }

        j = cr.End.J + 1;
        k = cr.Start.K - 1;
        u = U(i, j, k);
        u1 = U(i, j - 1, k);
        u2 = U(i, j, k + 1);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }

        j = cr.End.J + 1;
        k = cr.End.K + 1;
        u = U(i, j, k);
        u1 = U(i, j - 1, k);
        u2 = U(i, j, k - 1);
        for (int l = 0; l < dof; ++l)
        {
            u[l] = 0.5 * (u1[l] + u2[l]);
        }
    }
}

Structured<double>&
Block::UStorage(int i)
{
    switch (i)
    {
    case 0:
        return mU;
    case 1:
        return mU2;
    case 2:
        return mU3;
    default:
        throw 666;
    }
}

Structured<double>&
Block::UTStorage(int i)
{
    switch (i)
    {
    case 0:
        return mUT;
    case 1:
        return mUT2;
    case 2:
        return mUT3;
    default:
        throw 666;
    }
}

void
Block::ShiftTime()
{
    double* u3 = mU3.Data;
    mU3.Data = mU2.Data;
    mU2.Data = u3;
    mU2.SetTo(mU);

    double* ut3 = mUT3.Data;
    mUT3.Data = mUT2.Data;
    mUT2.Data = ut3;
    mUT2.SetTo(mUT);
}

