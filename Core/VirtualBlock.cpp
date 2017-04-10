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

#include "Roster.h"
#include "VirtualBlock.h"
#include "RigidBodyMotion.h"
#include <cstdlib>

VirtualBlock*
VirtualBlock::New(int id, int rank, const IndexRange& meshRange)
{
    VirtualBlock* block = new VirtualBlock(id, meshRange);
    Roster::GetInstance()->RegisterBlock(rank, block);

    return block;
}

VirtualBlock::VirtualBlock(int id, const IndexRange& meshRange)
:   mID(id),
    mMeshRange(meshRange),
    mCellRange(meshRange.Start + IndexIJK(1), meshRange.End),
    mRigidBodyMotion(NULL),
    mIsPeriodic(false), mRotPeriodicity(0.0, 0.0, 0.0)
{
}

VirtualBlock::~VirtualBlock()
{
    delete mRigidBodyMotion;
}

void
VirtualBlock::SetRigidBodyMotion(const RigidBodyMotion& motion)
{
    mRigidBodyMotion = motion.Clone();
}

bool
VirtualBlock::IsRotating() const
{
    if (mRigidBodyMotion != NULL)
    {
        RotationalMotion* rotmot = dynamic_cast<RotationalMotion*>(mRigidBodyMotion);
        if (rotmot != NULL)
        {
            return true;
        }
    }
    return false;
}

void
VirtualBlock::SetPeriodicity(const Vector3& rot)
{
    mRotPeriodicity = rot;
    mIsPeriodic = rot != Vector3(0.0, 0.0, 0.0);
}

