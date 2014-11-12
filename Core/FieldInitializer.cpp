// $Id: FieldInitializer.cpp 291 2013-07-25 10:25:11Z kato $

#include "FieldInitializer.h"
#include "Physics.h"
#include <string>
#include <cassert>

IsentropicExpansionInitializer::IsentropicExpansionInitializer(const FlowModel& model, const char* dir)
:   mModel(model), mDirection(X)
{
    std::string d(dir);
    if (d == "X")
        mDirection = X;
    else if (d == "Y")
        mDirection = Y;
    else if (d == "Z")
        mDirection = Z;
    else
        throw 666;
}

void
IsentropicExpansionInitializer::SetTotalQuantities(double P0, double T0)
{
    const Physics* PHYS = Physics::GetInstance();
    double pRef = PHYS->PRef();
    double TRef = PHYS->TRef();
    mP0 = P0 / pRef;
    mT0 = T0 / TRef;
}

void
IsentropicExpansionInitializer::AddPressureSpec(double pos, double pres)
{
    double pRef = Physics::GetInstance()->PRef();
    mSpecs.push_back(PressureSpec(pos, pres / pRef));
}

void
IsentropicExpansionInitializer::Initialize(Structured<double>& U, const Block& block) const
{
    IndexRange cr = block.CellRange();
    const Structured<double>& XYZ = block.XYZ();
    double gamma = Physics::GetInstance()->Gamma();
    double gm1 = gamma - 1.0;
    Vector3 nvel(0.0, 0.0, 0.0);
    nvel[mDirection] = 1.0;

    for (IndexIterator itor(cr); !itor.IsEnd(); itor.Advance())
    {
        IndexIJK ijk = itor.Index();
        double* u = U(ijk);
        double* xyz = XYZ(ijk);

        double phi, M, p, T, rho, c, v, rhoe, ke, rhoet;
        p = StaticPressureAt(xyz);

        phi = std::pow(mP0 / p, gm1 / gamma); // 1.0 + 0.5 * (gamma - 1.0) * M**2
        M = std::sqrt(2.0 * (phi - 1.0) / gm1);

        T = phi * mT0;
        rho = p / T;
        c = std::sqrt(gamma * p / rho);
        v = M * c;
        rhoe = p / gm1;
        ke = 0.5 * rho * v * v;
        rhoet = rhoe + ke;

        double UG[5];
        UG[0] = rho;
        UG[1] = rho * v * nvel.X();
        UG[2] = rho * v * nvel.Y();
        UG[3] = rho * v * nvel.Z();
        UG[4] = rhoet;

        mModel.FromGlobalToLocal(u, UG, block, ijk);

#if 1
        // Re-orient the velocity to align the specified direction *in the local frame*.
        v = std::sqrt(u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0];
        u[1] = u[0] * v * nvel.X();
        u[2] = u[0] * v * nvel.Y();
        u[3] = u[0] * v * nvel.Z();
#endif
    }
}

double
IsentropicExpansionInitializer::StaticPressureAt(double* xyz) const
{
    double pos = xyz[mDirection];
    double posMin = mSpecs.front().Position;
    double posMax = mSpecs.back().Position;
    pos = std::max(posMin, std::min(pos, posMax));

    size_t interval = 0;
    for (size_t i = 0; i < mSpecs.size() - 1; ++i)
    {
        if (mSpecs[i].Position <= pos && pos <= mSpecs[i + 1].Position)
        {
            interval = i;
            break;
        }
    }

    double pos1 = mSpecs[interval].Position;
    double pos2 = mSpecs[interval + 1].Position;
    double p1 = mSpecs[interval].Pressure;
    double p2 = mSpecs[interval + 1].Pressure;

    double xi = (pos - pos1) / (pos2 - pos1);
    return (1.0 - xi) * p1 + xi * p2;
}

CentrifugalInitializer::CentrifugalInitializer(const FlowModel& model, const char* axis)
:   mModel(model)
{
    std::string ax(axis);
    if (ax == "X")
        mAxis = X;
    else if (ax == "Y")
        mAxis = Y;
    else if (ax == "Z")
        mAxis = Z;
    else
        assert(false);
}

void
CentrifugalInitializer::SetTotalQuantities(double P0, double T0)
{
    const Physics* PHYS = Physics::GetInstance();
    double pRef = PHYS->PRef();
    double TRef = PHYS->TRef();
    mP0 = P0 / pRef;
    mT0 = T0 / TRef;
}

void
CentrifugalInitializer::AddSpec(double axialPos, double radialPos, const Vector3& velAR, double p)
{
    mSpecs.push_back(Spec(axialPos, radialPos, velAR, p));
}

void
CentrifugalInitializer::Initialize(Structured<double>& U, const Block& block) const
{
    IndexRange cr = block.CellRange();
    const Structured<double>& XYZ = block.XYZ();
    double gamma = Physics::GetInstance()->Gamma();
    double gm1 = gamma - 1.0;

    for (IndexIterator itor(cr); !itor.IsEnd(); itor.Advance())
    {
        IndexIJK ijk = itor.Index();
        double* u = U(ijk);
        double* xyz = XYZ(ijk);

        double phi, M, p, T, rho, c, v, rhoe, ke, rhoet;
        Vector3 nvel;
        GetSpecAt(xyz, nvel, p);

        phi = std::pow(mP0 / p, gm1 / gamma); // 1.0 + 0.5 * (gamma - 1.0) * M**2
        M = std::sqrt(2.0 * (phi - 1.0) / gm1);

        T = phi * mT0;
        rho = p / T;
        c = std::sqrt(gamma * p / rho);
        v = M * c;
        rhoe = p / gm1;
        ke = 0.5 * rho * v * v;
        rhoet = rhoe + ke;

        double UG[5];
        UG[0] = rho;
        UG[1] = rho * v * nvel.X();
        UG[2] = rho * v * nvel.Y();
        UG[3] = rho * v * nvel.Z();
        UG[4] = rhoet;

        mModel.FromGlobalToLocal(u, UG, block, ijk);

#if 0
        // Re-orient the velocity to align the specified direction *in the local frame*.
        v = std::sqrt(u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0];
        u[1] = u[0] * v * nvel.X();
        u[2] = u[0] * v * nvel.Y();
        u[3] = u[0] * v * nvel.Z();
#endif
    }
}

void
CentrifugalInitializer::GetSpecAt(double* xyz, Vector3& vel, double& pressure) const
{
    const double eps = 1.0e-5;
    double a, r;
    GetPosition(a, r, xyz);
    Vector3 p(a, r, 0.0);

    bool first = true;
    size_t imatch;
    double mindistSq, xi;
    for (size_t i = 0; i < mSpecs.size() - 1; ++i)
    {
        const Vector3& p1 = mSpecs[i].AR;
        const Vector3& p2 = mSpecs[i + 1].AR;
        Vector3 v12(p1, p2); // p2 - p1
        Vector3 v1p(p1, p); // p - p1
        double alpha = dot_product(v1p, v12) / v12.MagSq();
        if (alpha >= -eps && alpha <= (1.0 + eps))
        {
            Vector3 pi = p1 + alpha * v12;
            Vector3 d(pi, p);
            double distSq = d.MagSq();
            if (first || distSq < mindistSq)
            {
                mindistSq = distSq;
                imatch = i;
                xi = alpha;
                first = false;
            } 
        }
    }
    assert(!first);

    vel = (1.0 - xi) * mSpecs[imatch].VelAR + xi * mSpecs[imatch + 1].VelAR;
    pressure = (1.0 - xi) * mSpecs[imatch].Pressure + xi * mSpecs[imatch + 1].Pressure;
}

void
CentrifugalInitializer::GetPosition(double& axial, double& radial, double* xyz) const
{
    if (mAxis == X)
    {
        axial = xyz[0];
        radial = std::sqrt(xyz[1] * xyz[1] + xyz[2] * xyz[2]);
    }
    else if (mAxis == Y)
    {
        axial = xyz[1];
        radial = std::sqrt(xyz[0] * xyz[0] + xyz[2] * xyz[2]);
    }
    else if (mAxis == Z)
    {
        axial = xyz[2];
        radial = std::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
    }
    else
    {
        assert(false);
    }
}

