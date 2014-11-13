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

