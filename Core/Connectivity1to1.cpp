// $Id: Connectivity1to1.cpp 306 2013-10-02 07:03:25Z kato $

#include <mpi.h>
#include "Connectivity1to1.h"
#include "IndexUtils.h"
#include <cassert>

IndexTransform::IndexTransform(int t[3], const IndexRange& selfRange, const IndexRange& donorRange)
{
    IndexRange selfRangeFixed(selfRange), donorRangeFixed(donorRange);
    Canonize(t, selfRangeFixed, donorRangeFixed);
    SetTransformMatrix(t);
    Begin1 = selfRangeFixed.Start;
    Begin2 = donorRangeFixed.Start;
}

void
IndexTransform::Canonize(int t[3], IndexRange& selfRange, IndexRange& donorRange)
{
    if (selfRange.IsCanonical())
    {
        return;
    }

    std::cout << "IndexTransform: self-range is not in canonical form, fixing it" << std::endl;
    if (selfRange.Start.I > selfRange.End.I)
    {
        int i = std::abs(t[0]) - 1;
        int tmp;
        tmp = donorRange.Start[i];
        donorRange.Start[i] = donorRange.End[i];
        donorRange.End[i] = tmp;
    }
    if (selfRange.Start.J > selfRange.End.J)
    {
        int i = std::abs(t[1]) - 1;
        int tmp;
        tmp = donorRange.Start[i];
        donorRange.Start[i] = donorRange.End[i];
        donorRange.End[i] = tmp;
    }
    if (selfRange.Start.K > selfRange.End.K)
    {
        int i = std::abs(t[2]) - 1;
        int tmp;
        tmp = donorRange.Start[i];
        donorRange.Start[i] = donorRange.End[i];
        donorRange.End[i] = tmp;
    }

    selfRange.Canonize();
}

IndexTransform
IndexTransform::CellIndexTransform() const
{
    IndexTransform ci(*this);
    IndexIJK dBegin1(0, 0, 0);
    int d[3][3];
    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < 3; ++i)
        {
            d[i][j] = T[i][j] < 0 ? 1 : 0;
        }
        ci.Begin2.I += d[0][j];
        ci.Begin2.J += d[1][j];
        ci.Begin2.K += d[2][j];
    }

    return ci;
}

std::ostream&
IndexTransform::Dump(std::ostream& o) const
{
    o << "Begin1 = " << Begin1 << ", Begin2 = " << Begin2 << std::endl;
    for (int i = 0; i < 3; ++i)
    {
        o << '[';
        for (int j = 0; j < 3; ++j)
        {
            o << T[i][j] << ' ';
        }
        o << ']' << std::endl;
    }
    return o;
}

std::ostream& operator << (std::ostream& o, const IndexTransform& t)
{
    return t.Dump(o);
}

Connectivity1to1::Connectivity1to1(
    const IndexRange& meshRange, Direction direction, const Block& block,
    const IndexTransform& meshTransformSelfToDonor, int donorBlockID,
    int tag
    )
:   Connectivity(meshRange, direction),
    mBlockID(block.ID()),
    mDonorBlockID(donorBlockID),
    mTag(tag),
    mMeshTransformSelfToDonor(meshTransformSelfToDonor),
    mMeshTransformDonorToSelf(meshTransformSelfToDonor.Inverse()),
    mCellTransformSelfToDonor(meshTransformSelfToDonor.CellIndexTransform()),
    mPeriodicity(NONE)
{
    mCellRangeToSend =
        IndexUtils::FromMeshRangeToCellRange(
            meshRange, direction, block.GhostLayers()
            );
    mCellRangeToRecv =
        IndexUtils::FromMeshRangeToCellRange(
            meshRange, IndexUtils::Opposite(direction), block.GhostLayers()
            );
    mDonorCellRangeToRecv = mCellTransformSelfToDonor.Apply(mCellRangeToRecv).Canonical();
}

Connectivity1to1::~Connectivity1to1()
{
}

void
Connectivity1to1::SetPeriodicity(Periodicity periodicity, const Vector3& rotCenter, const Vector3& rotAngle)
{
    mPeriodicity = periodicity;
    mRotationCenter = rotCenter;
    mRotationAngle = rotAngle;
}

#include "Communicator.h"

void
Connectivity1to1::ApplyPeriodicity(Structured<double>& U) const // FIXME: what about turbulent quantities?
{
    if (!IsPeriodic())
    {
        return;
    }

    if (U.DOF() < 5)
    {
        // FIXME: this is very ugly.
        // if U's dof is less than 5, it probably doesn't contain conservative flow variables,
        // most likely turbulent quantities, so we skip this altogether.
        return;
    }

    Vector3 axis = mRotationAngle.Normalized();
    double ax = axis.X(), ay = axis.Y(), az = axis.Z();
    double theta = -mRotationAngle.Mag(); // The angle is from the self patch to the donor, so we need to negate this to rotate the received momentum into the self patch.
    double c = std::cos(theta), s = std::sin(theta);

    double R[3][3];
    R[0][0] = ax * ax + (1.0 - ax * ax) * c;
    R[0][1] = ax * ay * (1.0 - c) - az * s;
    R[0][2] = ax * az * (1.0 - c) + ay * s;
    R[1][0] = ax * ay * (1.0 - c) + az * s;
    R[1][1] = ay * ay + (1.0 - ay * ay) * c;
    R[1][2] = ay * az * (1.0 - c) - ax * s;
    R[2][0] = ax * az * (1.0 - c) - ay * s;
    R[2][1] = ay * az * (1.0 - c) + ax * s;
    R[2][2] = az * az + (1.0 - az * az) * c;

    int offset = 1;
    IndexRange cr = CellRangeToRecv();
    for (int k = cr.Start.K; k <= cr.End.K; ++k)
    {
        for (int j = cr.Start.J; j <= cr.End.J; ++j)
        {
            for (int i = cr.Start.I; i <= cr.End.I; ++i)
            {
                double* uu = U(i, j, k);
                Vector3 rhoV(&uu[offset]);
                rhoV.Apply(R);
                uu[offset    ] = rhoV.X();
                uu[offset + 1] = rhoV.Y();
                uu[offset + 2] = rhoV.Z();
            }
        }
    }
}

