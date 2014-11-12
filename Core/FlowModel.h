// $Id: FlowModel.h 284 2013-06-14 03:22:13Z kato $
#ifndef INCLUDED_FLOW_MODEL_H__
#define INCLUDED_FLOW_MODEL_H__

#include "ResidualEvaluator.h"
#include "Model.h"
#include "Structured.h"

//typedef ResidualEvaluator<5, Reconstructor<5, MinModLimiter, PrimitiveVariableCodec> > ResidualEval;

class Block;

class FlowModel : public Model
{
public:
    static void Initialize();

    FlowModel();

    int DOF() const { return 5; }

    Structured<double>& GetU(Block& block) const { return block.U(); }
    const Structured<double>& GetU(const Block& block) const { return block.U(); }

    void ApplyBCs(Block& block) const;

    double ComputeTimeStep(
        const Block& block,
        const Structured<double>& U,
        Structured<double>& DT,
        double cfl,
        bool localTimeStepping
        );

    // For matrix-free LUSGS
    // HH = 0.5 * (HL + HR) - 0.5 * |lambda| * (UR - UL)
    // AdU = (H(U + eps * dU) - H(U)) / eps
    // ScalarCoeff returns |lambda|
    // JacobianDU returns AdU

    double ScalarCoeff(
        const Block& block, const Structured<double>& U,
        const IndexIJK& Ii, const IndexIJK& Ij,
        double* Sij, double SijSign, double SijAbs,
        const Structured<double>& radius
        ) const;
    double ScalarCoeff(
        double* Ui, double* Uj, double* Sij, double SijSign, double SijAbs, double* Ri, double* Rj
        ) const;

    void JacobianDU(
        double* JacDU, const Block& block, const IndexIJK& i,
        double* U, double* dU, double* Sij, double SijSign, double SijAbs, double* Radius
        ) const;

#define SCALAR_DIAG 0
#if SCALAR_DIAG
    double DiagonalAddition(const Block& block, const Structured<double>& U, const IndexIJK& Ii) const { return 0.0; }
#else
    void DiagonalAddition(double* D, const Block& block, const Structured<double>& U, const IndexIJK& Ii) const { return; }
#endif

    // U + dU (but limits overshoots/undershoots)
    void UplusDU(Structured<double>& U, const Structured<double>& dU, const IndexRange& cellRange, const Block& block) const;

    // Coordinate frame transform
    virtual void FromGlobalToLocal(double* ULocal, double* UGlobal, const Block& block, const IndexIJK& ijk) const;
    virtual void FromLocalToGlobal(double* UGlobal, double* ULocal, const Block& block, const IndexIJK& ijk) const;

protected:

private:
    double Gamma;

    double mHartenEps; // FIXME
};

#endif // INCLUDED_FLOW_MODEL_H__

