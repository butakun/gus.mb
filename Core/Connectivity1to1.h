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
// $Id: Connectivity1to1.h 277 2013-06-04 01:58:51Z kato $
#ifndef INCLUDED_CONNECTIVITY_1_TO_1_H__
#define INCLUDED_CONNECTIVITY_1_TO_1_H__

#include "Connectivity.h"
#include "Vector3.h"
#include <cstdlib>
#include <iostream>

class IndexTransform
{
public:
    IndexTransform()
    {
    }

    IndexTransform(int t[3], const IndexIJK& begin1, const IndexIJK& begin2)
    : Begin1(begin1), Begin2(begin2)
    {
        SetTransformMatrix(t);
    }

    IndexTransform(int t[3], const IndexRange& selfRange, const IndexRange& donorRange);

    IndexTransform(const IndexTransform& rhs)
    : Begin1(rhs.Begin1), Begin2(rhs.Begin2)
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                T[i][j] = rhs.T[i][j];
    }

    IndexIJK Apply(const IndexIJK& i1) const
    {
        IndexIJK di1 = i1 - Begin1;
        return IndexIJK(
            T[0][0] * di1[0] + T[0][1] * di1[1] + T[0][2] * di1[2] + Begin2[0],
            T[1][0] * di1[0] + T[1][1] * di1[1] + T[1][2] * di1[2] + Begin2[1],
            T[2][0] * di1[0] + T[2][1] * di1[1] + T[2][2] * di1[2] + Begin2[2]
            );
    }

    IndexRange Apply(const IndexRange& range) const
    {
        return IndexRange(Apply(range.Start), Apply(range.End));
    }

    IndexTransform Inverse() const
    {
        IndexTransform Inv;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Inv.T[i][j] = T[j][i];
        Inv.Begin1 = Begin2;
        Inv.Begin2 = Begin1;
        return Inv;
    }

    IndexTransform CellIndexTransform() const;

    std::ostream& Dump(std::ostream& o) const;

    static void Canonize(int t[3], IndexRange& selfRange, IndexRange& donorRange);

protected:
    void SetTransformMatrix(int t[3])
    {
        T[0][0] = sign(t[0]) * del(t[0], 1); T[0][1] = sign(t[1]) * del(t[1], 1); T[0][2] = sign(t[2]) * del(t[2], 1);
        T[1][0] = sign(t[0]) * del(t[0], 2); T[1][1] = sign(t[1]) * del(t[1], 2); T[1][2] = sign(t[2]) * del(t[2], 2);
        T[2][0] = sign(t[0]) * del(t[0], 3); T[2][1] = sign(t[1]) * del(t[1], 3); T[2][2] = sign(t[2]) * del(t[2], 3);
    }

private:
    static int sign(int a) { return a >= 0 ? 1 : -1; }
    static int del(int a, int i) { return std::abs(a) == i ? 1 : 0; }

    IndexIJK Begin1, Begin2;
    int T[3][3];
};

std::ostream& operator << (std::ostream& o, const IndexTransform& t);

class Connectivity1to1 : public Connectivity
{
public:
    Connectivity1to1(
        const IndexRange& meshRange, Direction direction, const Block& block,
        const IndexTransform& meshTransformSelfToDonor, int donorBlockID,
        int tag
        );
    virtual ~Connectivity1to1();

    int BlockID() const { return mBlockID; }
    int DonorBlockID() const { return mDonorBlockID; }
    int Tag() const { return mTag; }

    IndexRange CellRangeToSend() const { return mCellRangeToSend; }
    IndexRange CellRangeToRecv() const { return mCellRangeToRecv; }
    IndexRange DonorCellRangeToRecv() const { return mDonorCellRangeToRecv; }

    IndexIJK DonorCellIndex(const IndexIJK& iSelf) const { return mCellTransformSelfToDonor.Apply(iSelf); }

    virtual void Apply(const Block& block, Structured<double>& U) {}
    virtual void ApplyTurb(const Block& block, Structured<double>& UT, const Structured<double>& U) {}

    // Periodicity
    typedef enum { NONE = 0, ROTATION, TRANSLATION } Periodicity;
    Periodicity GetPeriodicity() const { return mPeriodicity; }
    void SetPeriodicity(Periodicity periodicity, const Vector3& rotCenter, const Vector3& rotAngle);
    bool IsPeriodic() const { return mPeriodicity != NONE; }
    void ApplyPeriodicity(Structured<double>& U) const;

protected:

private:
    int mBlockID;
    int mDonorBlockID;
    int mTag;
    IndexRange mCellRangeToSend;
    IndexRange mCellRangeToRecv;
    IndexRange mDonorCellRangeToRecv;
    IndexTransform mMeshTransformSelfToDonor;
    IndexTransform mMeshTransformDonorToSelf;
    IndexTransform mCellTransformSelfToDonor;

    Periodicity mPeriodicity;
    // Rotational periodicity
    Vector3 mRotationCenter, mRotationAngle;
    // Translational periodicity
    Vector3 mTranslation;
};

#endif // INCLUDED_CONNECTIVITY_1_TO_1_H__

