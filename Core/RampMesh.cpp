// $Id: RampMesh.cpp 36 2010-06-07 00:30:59Z kato $

#include "RampMesh.h"
#include "Vector3.h"
#include <cmath>
#include <cassert>

void
RampMesh::GenerateMesh(Structured<double>& xyz, double LX, double LY, double LZ, double LX1, double delta1) const
{
    assert(xyz.Data != NULL);

    IndexRange range = xyz.GetRange();

    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        double zeta = double(k - range.Start.K) / double(range.End.K - range.Start.K);
        double z = LZ * zeta;
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            double eta = double(j - range.Start.J) / double(range.End.J - range.Start.J);
            eta = eta * eta;
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                double xi = double(i - range.Start.I) / double(range.End.I - range.Start.I);
                double x = LX * xi;
                double y, y0, y1;
                y0 = 0.0;
                if (x > LX1)
                {
                    y0 = (x - LX1) * std::atan(delta1 * M_PI / 180.0);
                }
                y1 = LY;
                y = (1.0 - eta) * y0 + eta * y1;
                xyz(i, j, k)[0] = x;
                xyz(i, j, k)[1] = y;
                xyz(i, j, k)[2] = z;
            }
        }
    }
}

