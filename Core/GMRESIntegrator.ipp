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
// $Id: GMRESIntegrator.ipp 116 2011-08-11 04:17:07Z kato $

template <class Model>
GMRESIntegrator<Model>::GMRESIntegrator(const Model& model, Block& block)
:   Integrator(model.DOF(), block, 1), // FIXME: eventually we need multiple steps to support multiple blocks / rank
    mModel(model),
    mRHS(mModel.DOF(), block.CellRange()),
    mAx(block, mRHS),
    mPreconditioner(block),
    mGMRES(25, mAx, mPreconditioner)
{
}

template <class Model>
GMRESIntegrator<Model>::~GMRESIntegrator()
{
    delete[] mRHS.Data;
}

