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
// $Id: BCFix.cpp 306 2013-10-02 07:03:25Z kato $

#include "BCFix.h"
#include "Physics.h"
#include "TurbulenceSpec.h"
#include <cassert>

BCFix::BCFix(
    const IndexRange& range, Direction direction,
    int ndof, double* ufix, int ndoft, const TurbulenceSpec* turbSpec
    )
:   BCPlanarLocal(range, direction),
    UFix(new double[ndof]),
    TurbFix(new double[ndoft])
{
    Physics* PHYS = Physics::GetInstance();
    double Uref[5];
    Uref[0] = PHYS->RhoRef();
    Uref[1] = PHYS->RhoRef() * PHYS->VRef();
    Uref[2] = PHYS->RhoRef() * PHYS->VRef();
    Uref[3] = PHYS->RhoRef() * PHYS->VRef();
    Uref[4] = PHYS->RhoRef() * PHYS->ERef();

    for (int i = 0; i < ndof; ++i)
    {
        UFix[i] = ufix[i] / Uref[i];
    }

    //double UTref[2];
    //UTref[0] = PHYS->RhoRef() * PHYS->TKERef(); // rho * tke
    //UTref[1] = PHYS->RhoRef() * PHYS->OmegaRef(); // rho * omega
    assert(ndoft == 2);

    // FIXME: how do we make this model-neutral?
    TurbFix[0] = turbSpec->TKE_Nondimensional();
    TurbFix[1] = turbSpec->Omega_Nondimensional();
}

BCFix::~BCFix()
{
    delete[] UFix;
    delete[] TurbFix;
}

void
BCFix::SetTo(int ndof, double* ufix, int ndoft, const TurbulenceSpec* turbSpec)
{
    Physics* PHYS = Physics::GetInstance();
    double Uref[5];
    Uref[0] = PHYS->RhoRef();
    Uref[1] = PHYS->RhoRef() * PHYS->VRef();
    Uref[2] = PHYS->RhoRef() * PHYS->VRef();
    Uref[3] = PHYS->RhoRef() * PHYS->VRef();
    Uref[4] = PHYS->RhoRef() * PHYS->ERef();

    for (int i = 0; i < ndof; ++i)
    {
        UFix[i] = ufix[i] / Uref[i];
    }

    // FIXME: how do we make this model-neutral?
    TurbFix[0] = turbSpec->TKE_Nondimensional();
    TurbFix[1] = turbSpec->Omega_Nondimensional();
}

void
BCFix::LocalFunc(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& U, const Block& block
    )
{
    IndexIJK ijk = iGhost;
    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* u = U(ijk);
        for (int l = 0; l < U.DOF(); ++l) // FIXME: what if U.NDOF is different from UFix dof?
        {
            u[l] = UFix[l];
        }
        ijk = ijk - deltaInterior;
    }
}

void
BCFix::LocalFuncTurb(
    const IndexIJK& iFace,
    const IndexIJK& iGhost, const IndexIJK& iInterior, const IndexIJK& deltaInterior,
    const Vector3& Sn, Structured<double>& UT, const Structured<double>& U, const Block& block
    )
{
    IndexIJK ijk = iGhost;
    for (int i = 0; i < block.GhostLayers(); ++i)
    {
        double* ut = UT(ijk);
        double* u = U(ijk);
        double rho = u[0];
        // FIXME: model-specific
        ut[0] = rho * TurbFix[0];
        ut[1] = rho * TurbFix[1];
        ijk = ijk - deltaInterior;
    }
}

