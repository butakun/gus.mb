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
// $Id: Test.cpp 116 2011-08-11 04:17:07Z kato $

#include "Physics.h"
#include "Block.h"
#include "RampMesh.h"
#include "ResidualEvaluator.h"
#include "SolverUtils.h"
#include "Vector3.h"
#include "FlowModel.h"
#include "SimpleExplicitIntegrator.h"
#include "LUSGSIntegrator.h"
#include <cstdlib>
#include <iostream>

int main()
{
    IndexRange meshRange(0, 0, 0, 10, 1, 1);
    std::cout << meshRange.Size() << std::endl;

    Block block(1, meshRange);
    IndexRange cellRange = block.CellRange();

    RampMesh mesher;
    mesher.GenerateMesh(block.XYZ(), 20.0, 10.0, 1.0, 10.0, 10.0);
    block.ComputeMetrics();

    double U0[5] = { 1.0, 2.0, 0.0, 0.0, 4.0 };
    block.U().SetTo(U0);

    double* v = block.U()(5, 0, 0);
    std::cout << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ", " << v[4] << std::endl;

    Structured<double> R(5, block.Dim());
    Structured<double> Lambda(1, block.Dim());
    SolverUtils::ComputeSpectralRadius(block, block.U(), Lambda);

    std::cout << "start" << std::endl;
    R.SetTo(0.0);
    ResidualEvaluator resEval;
    resEval.EvaluateResidual(block, block.U(), R);

    for (int i = cellRange.Start.I; i <= cellRange.End.I; ++i)
    {
        int j = cellRange.Start.J;
        int k = cellRange.Start.K;
        double* r = R(i, j, k);
        Vector3 sxi(block.Sxi()(i, j, k));
        Vector3 seta(block.Seta()(i, j, k));
        Vector3 szeta(block.Szeta()(i, j, k));
        std::cout << i << ": " << r[0] << ", " << r[1] << ", " << r[2] << ", " << r[3] << ", " << r[4] << ", lambda = " << *Lambda(i, j, k) << std::endl;
    }
    std::cout << "end" << std::endl;

#if 0
    SimpleExplicitIntegrator integrator(block);
    Structured<double>& U = block.U();
    integrator.Integrate(0, U);
#else
    FlowModel model;
    LUSGSIntegrator<FlowModel> integrator(model, block, block.U());
    integrator.Integrate(0, block.U());
    integrator.Integrate(1, block.U());
#endif

    return EXIT_SUCCESS;
}

