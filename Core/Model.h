// $Id: Model.h 253 2012-07-20 10:03:41Z kato $
#ifndef INCLUDED_MODEL_H__
#define INCLUDED_MODEL_H__

#include "Block.h"

class Model
{
public:
    virtual ~Model() {}

    // Coordinate frame transform
    virtual void FromGlobalToLocal(double* ULocal, double* UGlobal, const Block& block, const IndexIJK& ijk) const = 0;
    virtual void FromLocalToGlobal(double* UGlobal, double* ULocal, const Block& block, const IndexIJK& ijk) const = 0;
};

#endif // INCLUDED_MODEL_H__

