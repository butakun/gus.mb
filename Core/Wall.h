// $Id: Wall.h 37 2010-06-11 16:00:18Z kato $
#ifndef INCLUDED_WALL_H__
#define INCLUDED_WALL_H__

#include "Structured.h"

class Wall
{
public:
    Wall(const IndexRange& meshRange)
    :   mXYZ(3, meshRange)
    {}

    ~Wall()
    {
        delete[] mXYZ.Data;
    }

    const Structured<double>& XYZ() const { return mXYZ; }

protected:

private:
    Structured<double> mXYZ;
};

#endif // INCLUDED_WALL_H__

