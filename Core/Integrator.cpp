// $Id: Integrator.cpp 60 2010-09-01 16:15:57Z kato $

#include "Integrator.h"

Integrator::Integrator(int dof, Block& block, int numSteps)
:   mBlock(block), mIntegrationSteps(numSteps),
    mDU(dof, block.CellRange()),
    mDT(1, block.CellRange())
{
}

Integrator::~Integrator()
{
    delete[] mDU.Data;
    delete[] mDT.Data;
}

