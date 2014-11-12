// $Id: BCFunctors.h 4 2010-03-06 14:10:00Z kato $
#ifndef INCLUDED_HIRO_BC_FUNCTORS_H__
#define INCLUDED_HIRO_BC_FUNCTORS_H__

#include "Block.h"
#include "Structured.h"

class BCFunctor
{
public:
    virtual ~BCFunctor() {}

    virtual void Apply(
        Structured<double> U,
        const Block& block
        ) = 0;

protected:

private:
};

class BCInviscidWall : public BCFunctor
{
public:
    BCInviscidWall(const IndexRange& range, int direction);
    virtual ~BCInviscidWall();

    virtual void Apply(
        Structured<double> U,
        const Block& block
        );

protected:

private:
    IndexRange mRange;
    int mDirection;
};

class BCExtrapolate : public BCFunctor
{
public:
    BCExtrapolate(int iin, int ighost);
    virtual ~BCExtrapolate();

    virtual void Apply(
        Structured<double> U,
        const Block& block
        );

protected:

private:
    int mIIn, mIGhost;
};

#endif // INCLUDED_HIRO_BC_FUNCTORS_H__

