// $Id: VirtualBlock.cpp 244 2012-06-01 15:39:32Z kato $

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

