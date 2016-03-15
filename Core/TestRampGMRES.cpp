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
// $Id$

#include "Block.h"
#include "Physics.h"
#include "Communicator.h"
#include "RampMesh.h"
#include "ResidualEvaluator.h"

class Operator
{
public:
    Operator(const Block& block) : mBlock(block) {}

protected:

private:
    Structured<double> U0;
    Structured<double> MuK0;
    Structured<double> R0, R;
};

class BlockVector
{
public:

protected:

private:
};

int main()
{
    Communicator::Initialize(&argc, &argv);

#if 0
    DomainInfo di;
    Domain domain(di);
#endif

    int imin = 0, jmin = 0, kmin = 0, imax = 100, jmax = 50, kmax = 1;
    IndexRange meshRange(imin, jmin, kmin, imax, jmax, kmax);
    Block block(1, meshRange);
    domain.RegisterLocalBlock(&block);

    RampMesh mesher;
    mesher.GenerateMesh(block.XYZ(), 20.0, 10.0, 1.0, 10.0, 10.0);
    block.ComputeMetrics();

    double M0 = 2.5, Rho0 = 1.2, T0 = 300.0, gamma, RGAS, C0, V0, Et;
    Physics::Initialize(1.4, Rho0, T0);
    gamma = Physics::GetInstance()->Gamma();
    RGAS = Physics::GetInstance()->RGAS();
    C0 = std::sqrt(gamma * Physics::GetInstance()->RGAS() * T0);
    V0 = M0 * C0;
    Et = RGAS / (gamma - 1.0) * T0 + 0.5 * V0 * V0;
    double rho0, v0, rhoet0;
    rho0 = Rho0 / Physics::GetInstance()->RhoRef();
    v0 = V0 / Physics::GetInstance()->VRef();
    rhoet0 = Rho0 * Et / (Physics::GetInstance()->RhoRef() * Physics::GetInstance()->ERef());
    std::cout << "rho0, v0, rhoet0 = " << rho0 << ", " << v0 << ", " << rhoet0 << std::endl;

    double U0[5] = { rho0, v0, 0.0, 0.0, rhoet0 };
    double UT0[2] = { 0.0, 1.0 };
    block.U().SetTo(U0);
    block.UT().SetTo(UT0);

}

