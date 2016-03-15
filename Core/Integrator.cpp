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

