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

