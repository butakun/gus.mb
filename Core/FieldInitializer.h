/*
    gus.mb, an open source flow solver.
    Copyright (C) 2016 Hiromasa Kato <hiromasa at gmail.com>

    This file is part of gus.mb.

    gus.mb is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    gus.mb is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
// $Id: FieldInitializer.h 291 2013-07-25 10:25:11Z kato $
#ifndef INCLUDED_FIELD_INITIALIZER_H__
#define INCLUDED_FIELD_INITIALIZER_H__

#include "FlowModel.h"
#include "Structured.h"
#include "Block.h"

class FieldInitializer
{
public:
    enum Direction {X = 0, Y = 1, Z = 2};

    virtual ~FieldInitializer() {}
    virtual void Initialize(Structured<double>& U, const Block& block) const = 0;

protected:

private:
};

class IsentropicExpansionInitializer : public FieldInitializer
{
public:
    IsentropicExpansionInitializer(const FlowModel& model, const char* dir);
    virtual ~IsentropicExpansionInitializer() {}

    virtual void Initialize(Structured<double>& U, const Block& block) const;

    void SetTotalQuantities(double P0, double T0); // dimensional values expected.
    void AddPressureSpec(double pos, double pres); // dimensional pressure expected.

protected:
    double StaticPressureAt(double* xyz) const;

private:
    class PressureSpec
    {
    public:
        double Position, Pressure; // Pressure is nondimensional.

        PressureSpec(double pos = 0.0, double pres = 0.0) : Position(pos), Pressure(pres) {}
    };

    FlowModel mModel;
    Direction mDirection;

    double mP0, mT0; // Pressure/temperature are both nondimensional.
    std::vector<PressureSpec> mSpecs;
};

class CentrifugalInitializer : public FieldInitializer
{
public:
    CentrifugalInitializer(const FlowModel& model, const char* axis);
    virtual ~CentrifugalInitializer() {}

    virtual void Initialize(Structured<double>& U, const Block& block) const;

    void SetTotalQuantities(double P0, double T0); // dimensional values expected.
    void AddSpec(double axialPos, double radialPos, const Vector3& velAR, double p);

protected:

private:
    class Spec
    {
    public:
        Vector3 AR;
        Vector3 VelAR;
        double Pressure;
        Spec(double axial, double radial, const Vector3& velAR, double p)
        :   AR(axial, radial, 0.0), VelAR(velAR), Pressure(p)
        {
            VelAR.Normalize();
        }
    };

    void GetPosition(double& axial, double& radial, double* xyz) const;
    void GetSpecAt(double* xyz, Vector3& vel, double& pressure) const;

    FlowModel mModel;
    Direction mAxis;
    double mP0, mT0;
    std::vector<Spec> mSpecs;
};

#endif // INCLUDED_FIELD_INITIALIZER_H__

