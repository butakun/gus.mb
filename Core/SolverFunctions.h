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
// $Id: SolverFunctions.h 4 2010-03-06 14:10:00Z kato $
#ifndef INCLUDED_HIRO_SOLVER_FUNCTIONS_H__
#define INCLUDED_HIRO_SOLVER_FUNCTIONS_H__

#include "Structured.h"
#include "Block.h"

void EvaluateResidual(Block& block, Structured<double> R, int order);

void EvaluateResidual(Block& block, Structured<double> U, Structured<double> R, int order);

void UpdateSolution(
    IndexRange cellRange,
    Structured<double> U,
    Structured<double> DT,
    Structured<double> Vol,
    Structured<double> RHS
    );

void EvaluateFlux(
    double* UL, double* UR,
    double* Sn,
    double gamma,
    double* H
    );

void EvaluateResidual(
    double* U,
    double* Sxi, double* Seta,
    int iStart, int iEnd,
    int jStart, int jEnd,
    int idim, int jdim,
    double* R
    );

void ComputeMetrics(
    double* XYZ,
    double* Sxi, double* Seta, double* Vol,
    int istart, int iend,
    int jstart, int jend,
    int idim, int jdim
    );

void ComputeSpectralRadius(
    double* U,
    double* Sxi, double* Seta, double* Vol,
    double gamma,
    int istart, int iend,
    int jstart, int jend,
    int idim, int jdim,
    double* Lambda
    );

#endif // INCLUDED_HIRO_SOLVER_FUNCTIONS_H__

