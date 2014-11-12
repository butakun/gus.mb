// $Id: Connectivity.h 44 2010-07-08 16:15:52Z kato $
#ifndef INCLUDED_CONNECTIVITY_H__
#define INCLUDED_CONNECTIVITY_H__

#include "BCPlanar.h"

class Connectivity : public BCPlanar
{
public:
    Connectivity(const IndexRange& meshRange, Direction direction) : BCPlanar(meshRange, direction) {}
    virtual ~Connectivity() {}

protected:

private:
};

#endif // INCLUDED_CONNECTIVITY_H__

