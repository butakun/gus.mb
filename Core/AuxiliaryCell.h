// $Id: AuxiliaryCell.h 60 2010-09-01 16:15:57Z kato $
#ifndef INCLUDED_AUXILIARY_CELL_H__
#define INCLUDED_AUXILIARY_CELL_H__

#include "Block.h"
#include "Structured.h"
#include "IndexUtils.h"
#include <cassert>

/// Adapter class representing a dual cell, used for the evaluation of gradient at flux face.
class AuxiliaryCell
{
public:
    // iU -> see the documentation (ViscousFlux.fig)
    AuxiliaryCell(
        Direction dir,
        const IndexIJK& iU_,
        const Structured<double>& Sxi_,
        const Structured<double>& Seta_,
        const Structured<double>& Szeta_,
        const Structured<double>& Vol_,
        const IndexRange& cellRange
        );

    AuxiliaryCell(Direction dir, const IndexIJK& iU_, const Block& block);

    void SetIndex(const IndexIJK& i) { Init(i); }
    const IndexIJK& IndexU1M() const { return iU1M; }
    const IndexIJK& IndexU1P() const { return iU1P; }
    void IndicesU2M(IndexIJK& i1, IndexIJK& i2, IndexIJK& i3, IndexIJK& i4) const;
    void IndicesU2P(IndexIJK& i1, IndexIJK& i2, IndexIJK& i3, IndexIJK& i4) const;
    void IndicesU3M(IndexIJK& i1, IndexIJK& i2, IndexIJK& i3, IndexIJK& i4) const;
    void IndicesU3P(IndexIJK& i1, IndexIJK& i2, IndexIJK& i3, IndexIJK& i4) const;

    const Vector3& Surface() const { return Sn; } // This is the cell face around which this auxiliary cell is defined.
    const Vector3& Surface1M() const { return S1M; }
    const Vector3& Surface1P() const { return S1P; }
    const Vector3& Surface2M() const { return S2M; }
    const Vector3& Surface2P() const { return S2P; }
    const Vector3& Surface3M() const { return S3M; }
    const Vector3& Surface3P() const { return S3P; }
    double Volume() const { return Volu; }

protected:

private:
    void Init(const IndexIJK& i);

    Direction Dir;
    IndexIJK iU;
    IndexRange CellRange;
    const Structured<double>& Sxi;
    const Structured<double>& Seta;
    const Structured<double>& Szeta;
    const Structured<double>& Vol;

    IndexIJK iUD1, iUD2, iUD3;
    IndexIJK iU1M, iU1P;

    Vector3 Sn;
    Vector3 S1M, S1P, S2M, S2P, S3M, S3P;
    double Volu;
};

inline
AuxiliaryCell::AuxiliaryCell(
    Direction dir,
    const IndexIJK& iU_,
    const Structured<double>& Sxi_,
    const Structured<double>& Seta_,
    const Structured<double>& Szeta_,
    const Structured<double>& Vol_,
    const IndexRange& cellRange
    )
:   Dir(dir), iU(iU_), CellRange(cellRange), Sxi(Sxi_), Seta(Seta_), Szeta(Szeta_), Vol(Vol_)
{
    Init(iU);
}

inline
AuxiliaryCell::AuxiliaryCell(Direction dir, const IndexIJK& iU_, const Block& block)
:   Dir(dir), iU(iU_), CellRange(block.CellRange()), Sxi(block.Sxi()), Seta(block.Seta()), Szeta(block.Szeta()), Vol(block.Vol())
{
    Init(iU);
}

inline
void
AuxiliaryCell::Init(const IndexIJK& i)
{
    iU = i;

    const Structured<double> *S1 = NULL, *S2 = NULL, *S3 = NULL;
    bool boundaryM = false, boundaryP = false;

    switch (Dir)
    {
    case I:
        iUD1 = IndexIJK(1, 0, 0);
        iUD2 = IndexIJK(0, 1, 0);
        iUD3 = IndexIJK(0, 0, 1);
        S1 = &Sxi;
        S2 = &Seta;
        S3 = &Szeta;
        boundaryM = iU.I == CellRange.Start.I - 1;
        boundaryP = iU.I == CellRange.End.I;
        break;
    case J:
        iUD1 = IndexIJK(0, 1, 0);
        iUD2 = IndexIJK(0, 0, 1);
        iUD3 = IndexIJK(1, 0, 0);
        S1 = &Seta;
        S2 = &Szeta;
        S3 = &Sxi;
        boundaryM = iU.J == CellRange.Start.J - 1;
        boundaryP = iU.J == CellRange.End.J;
        break;
    case K:
        iUD1 = IndexIJK(0, 0, 1);
        iUD2 = IndexIJK(1, 0, 0);
        iUD3 = IndexIJK(0, 1, 0);
        S1 = &Szeta;
        S2 = &Sxi;
        S3 = &Seta;
        boundaryM = iU.K == CellRange.Start.K - 1;
        boundaryP = iU.K == CellRange.End.K;
        break;
    default:
        assert(false);
    }

    iU1M = iU;
    iU1P = iU + iUD1;

    // Surface normals
    IndexIJK iM1M, iM1, iM1P;
    iM1M = iU - iUD1;
    iM1  = iU;
    iM1P = iU + iUD1;

    IndexIJK iM211, iM212, iM221, iM222, iM311, iM312, iM321, iM322;
    iM211 = iU - iUD2;
    iM212 = iU + iUD1 - iUD2;
    iM221 = iU;
    iM222 = iU + iUD1;
    iM311 = iU - iUD3;
    iM312 = iU + iUD1 - iUD3;
    iM321 = iU;
    iM322 = iU + iUD1;

    if (boundaryM)
    {
        iM1M  = iM1;
        iM211 = iM212;
        iM221 = iM222;
        iM311 = iM312;
        iM321 = iM322;
    }
    else if (boundaryP)
    {
        iM1P  = iM1;
        iM212 = iM211;
        iM222 = iM221;
        iM312 = iM311;
        iM322 = iM321;
    }

    Sn = Vector3((*S1)(iM1));
    S1M = 0.5 * (Vector3((*S1)(iM1)) + Vector3((*S1)(iM1M)));
    S1P = 0.5 * (Vector3((*S1)(iM1)) + Vector3((*S1)(iM1P)));
    S2M = 0.5 * (Vector3((*S2)(iM211)) + Vector3((*S2)(iM212)));
    S2P = 0.5 * (Vector3((*S2)(iM221)) + Vector3((*S2)(iM222)));
    S3M = 0.5 * (Vector3((*S3)(iM311)) + Vector3((*S3)(iM312)));
    S3P = 0.5 * (Vector3((*S3)(iM321)) + Vector3((*S3)(iM322)));

    Volu = 0.5 * (*Vol(iU) + *Vol(iU + iUD1));
}

inline
void
AuxiliaryCell::IndicesU2M(IndexIJK& i1, IndexIJK& i2, IndexIJK& i3, IndexIJK& i4) const
{
    i1 = iU1M;
    i2 = iU1P;
    i3 = iU - iUD2;
    i4 = iU + iUD1 - iUD2;
}

inline
void
AuxiliaryCell::IndicesU2P(IndexIJK& i1, IndexIJK& i2, IndexIJK& i3, IndexIJK& i4) const
{
    i1 = iU1M;
    i2 = iU1P;
    i3 = iU + iUD2;
    i4 = iU + iUD1 + iUD2;
}

inline
void
AuxiliaryCell::IndicesU3M(IndexIJK& i1, IndexIJK& i2, IndexIJK& i3, IndexIJK& i4) const
{
    i1 = iU1M;
    i2 = iU1P;
    i3 = iU - iUD3;
    i4 = iU + iUD1 - iUD3;
}

inline
void
AuxiliaryCell::IndicesU3P(IndexIJK& i1, IndexIJK& i2, IndexIJK& i3, IndexIJK& i4) const
{
    i1 = iU1M;
    i2 = iU1P;
    i3 = iU + iUD3;
    i4 = iU + iUD1 + iUD3;
}

#endif // INCLUDED_AUXILIARY_CELL_H__

