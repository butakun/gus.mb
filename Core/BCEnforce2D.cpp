// $Id: BCEnforce2D.cpp 175 2012-01-04 06:01:45Z kato $

#include "BCEnforce2D.h"

BCEnforce2D::BCEnforce2D(Type type)
:   BC(IndexRange(-1, -1, -1, -1, -1, -1)), mType(type)
{
}

BCEnforce2D::~BCEnforce2D()
{
}

void
BCEnforce2D::Apply(const Block& block, Structured<double>& U)
{
    IndexRange r = U.GetRange();
    for (IndexIterator it(r); !it.IsEnd(); it.Advance())
    {
        IndexIJK ijk = it.Index();

        // FIXME: only TWOD_Z is implemented.
        double* UU = U(ijk);
        Vector3 e(UU[1], UU[2], 0.0);
        e.Normalize();
        double rhoV = std::sqrt(UU[1] * UU[1] + UU[2] * UU[2] + UU[3] * UU[3]);
        UU[1] = rhoV * e.X();
        UU[2] = rhoV * e.Y();
        UU[3] = rhoV * e.Z();
    }
}

void
BCEnforce2D::ApplyTurb(const Block& block, Structured<double>& UT)
{
}

