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
// $Id: VirtualBlock.h 258 2012-12-26 07:36:54Z kato $
#ifndef INCLUDED_VIRTUAL_BLOCK_H__
#define INCLUDED_VIRTUAL_BLOCK_H__

#include "Structured.h"
#include "Vector3.h"

class RigidBodyMotion;

class VirtualBlock
{
public:
    static VirtualBlock* New(int id, int rank, const IndexRange& meshRange);

    virtual ~VirtualBlock();

    int ID() const { return mID; }

    const IndexRange& MeshRange() const { return mMeshRange; }
    const IndexRange& CellRange() const { return mCellRange; }

    void SetRigidBodyMotion(const RigidBodyMotion& motion);
    RigidBodyMotion* GetRigidBodyMotion() const { return mRigidBodyMotion; }
    bool IsStationary() const { return mRigidBodyMotion == NULL; }
    bool IsRotating() const;

    bool IsPeriodic() const { return mIsPeriodic; }
    void SetPeriodicity(const Vector3& rot);
    const Vector3& Periodicity() const { return mRotPeriodicity; }

protected:
    VirtualBlock(int id, const IndexRange& meshRange);

private:
    int mID;
    IndexRange mMeshRange;
    IndexRange mCellRange;  // the index range covering the interior cells, i.e. not ghost cells.
    RigidBodyMotion* mRigidBodyMotion; // NULL if this block is stationary.
    bool mIsPeriodic;
    Vector3 mRotPeriodicity; // FIXME: what about translational periodicity?
};

#endif // INCLUDED_VIRTUAL_BLOCK_H__

