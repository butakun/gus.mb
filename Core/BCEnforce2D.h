// $Id: BCEnforce2D.h 175 2012-01-04 06:01:45Z kato $
#ifndef INCLUDED_ENFORCE_2D_H__
#define INCLUDED_ENFORCE_2D_H__

#include "BC.h"

class BCEnforce2D : public BC
{
public:
    enum Type { TWOD_X, TWOD_Y, TWOD_Z, AXISYM_X, AXISYM_Y, AXISYM_Z };

    BCEnforce2D(Type type);
    virtual ~BCEnforce2D();

    virtual void Apply(const Block& block, Structured<double>& U);
    virtual void ApplyTurb(const Block& block, Structured<double>& UT);

protected:

private:
    Type mType;
};

#endif // INCLUDED_ENFORCE_2D_H__

