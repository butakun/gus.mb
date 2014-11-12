// $Id: Structured.cpp 14 2010-04-14 10:56:36Z kato $

#include "Structured.h"

IndexRange
IndexRange::Canonical() const
{
    int imin, jmin, kmin, imax, jmax, kmax;
    imin = std::min(Start.I, End.I);
    imax = std::max(Start.I, End.I);
    jmin = std::min(Start.J, End.J);
    jmax = std::max(Start.J, End.J);
    kmin = std::min(Start.K, End.K);
    kmax = std::max(Start.K, End.K);

    return IndexRange(imin, jmin, kmin, imax, jmax, kmax);
}

