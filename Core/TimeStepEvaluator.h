// $Id: TimeStepEvaluator.h 255 2012-12-05 07:47:00Z kato $
#ifndef INCLUDED_TIME_STEP_EVALUATOR_H__
#define INCLUDED_TIME_STEP_EVALUATOR_H__

#include "Integrator.h"
#include "Communicator.h"
#include <vector>
#include <mpi.h>

template <class Model>
class TimeStepEvaluator
{
public:
    TimeStepEvaluator(Model model, double cfl, bool localTimeStepping)
    : mModel(model), mCFL(cfl), mLocalTimeStepping(localTimeStepping)
    {
    }

    void Evaluate(Integrators& integrators);

    void SetCFL(double cfl) { mCFL = cfl; }

protected:

private:
    Model mModel;
    double mCFL;
    bool mLocalTimeStepping;
};

template <class Model>
void
TimeStepEvaluator<Model>::Evaluate(Integrators& integrators)
{
    Communicator* COMM = Communicator::GetInstance();

    double dtmin;
    bool first = true;
    for (Integrators::iterator i = integrators.begin();
        i != integrators.end(); ++i)
    {
        Integrator* integrator = *i;
        const Block& block = integrator->GetBlock();

        double dtminLocal;
        dtminLocal =
            mModel.ComputeTimeStep(
                block,
                mModel.GetU(block),
                integrator->DT(),
                mCFL,
                mLocalTimeStepping
                );
        if (first)
        {
            dtmin = dtminLocal;
            first = false;
        }
        else
        {
            dtmin = std::min(dtminLocal, dtmin);
        }
    }

    // Gather from all other processes
    int err;
    double* dtmins = new double[COMM->Size()];
    err = MPI_Allgather(&dtmin, 1, MPI_DOUBLE, dtmins, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
    for (int i = 0; i < COMM->Size(); ++i)
    {
        dtmin = std::min(dtmin, dtmins[i]);
    }
    delete[] dtmins;

    // At this point we now have the global minimal value of dt.
    // We then limit the dt at most C * dtmin in all the blocks;

    double C, dtmax;
    if (mLocalTimeStepping)
    {
        C = 5.0;
    }
    else
    {
        C = 1.0;
    }
    dtmax = C * dtmin;

    for (Integrators::iterator i = integrators.begin();
        i != integrators.end(); ++i)
    {
        Integrator* integrator = *i;
        Structured<double>& DT = integrator->DT();
        const Block& block = integrator->GetBlock();
        IndexRange cr = block.CellRange();
        for (int k = cr.Start.K; k <= cr.End.K; ++k)
        {
            for (int j = cr.Start.J; j <= cr.End.J; ++j)
            {
                for (int i = cr.Start.I; i <= cr.End.I; ++i)
                {
                    double dt = *DT(i, j, k);
                    *DT(i, j, k) = std::min(dt, dtmax);
                }
            }
        }
    }
}

#endif // INCLUDED_TIME_STEP_EVALUATOR_H__

