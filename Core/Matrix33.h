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
// $Id: Matrix33.h 312 2013-12-13 08:50:20Z kato $
#ifndef INCLUDED_MATRIX_33_H__
#define INCLUDED_MATRIX_33_H__

#include "Vector3.h"
#include <iostream>

class Matrix33
{
public:
    double M[3][3];

    Matrix33() {};
    Matrix33(const Matrix33& m);

    static Matrix33 RotationMatrix(const Vector3& rotation);

    double operator () (int i, int j) const { return M[i][j]; }
    double& operator () (int i, int j) { return M[i][j]; }

    void Multiply(Vector3& v) const; // v = M * v
    void Multiply(double* v) const; // v = M * v

    Matrix33 Inverse() const;

protected:

private:
};

inline
Matrix33::Matrix33(const Matrix33& m)
{
    M[0][0] = m(0, 0); M[0][1] = m(0, 1); M[0][2] = m(0, 2);
    M[1][0] = m(1, 0); M[1][1] = m(1, 1); M[1][2] = m(1, 2);
    M[2][0] = m(2, 0); M[2][1] = m(2, 1); M[2][2] = m(2, 2);
}

inline
Matrix33
Matrix33::RotationMatrix(const Vector3& rot)
{
    Vector3 u(0.0, 0.0, 0.0);
    if (rot != Vector3(0.0, 0.0, 0.0))
    {
        u = rot.Normalized();
    }
    double theta = rot.Mag();
    double cos = std::cos(theta), sin = std::sin(theta), onemcos = 1.0 - cos;

    Matrix33 R;

    R(0, 0) = cos + u.X() * u.X() * onemcos;
    R(1, 0) = u.Y() * u.X() * onemcos + u.Z() * sin;
    R(2, 0) = u.Z() * u.X() * onemcos - u.Y() * sin;

    R(0, 1) = u.X() * u.Y() * onemcos - u.Z() * sin;
    R(1, 1) = cos + u.Y() * u.Y() * onemcos;
    R(2, 1) = u.Z() * u.Y() * onemcos + u.X() * sin;

    R(0, 2) = u.X() * u.Z() * onemcos + u.Y() * sin;
    R(1, 2) = u.Y() * u.Z() * onemcos - u.X() * sin;
    R(2, 2) = cos + u.Z() * u.Z() * onemcos;

    return R;
}

inline
void
Matrix33::Multiply(Vector3& v) const
{
    double r[3];
    r[0] = M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2];
    r[1] = M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2];
    r[2] = M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2];
    v[0] = r[0];
    v[1] = r[1];
    v[2] = r[2];
}

inline
void
Matrix33::Multiply(double* v) const
{
    double r[3];
    r[0] = M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2];
    r[1] = M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2];
    r[2] = M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2];
    v[0] = r[0];
    v[1] = r[1];
    v[2] = r[2];
}

inline
Matrix33
Matrix33::Inverse() const
{
    Matrix33 Inv;
    double d11, d12, d13;
    double d21, d22, d23;
    double d31, d32, d33;
    double deti;

    d11 = M[1][1] * M[2][2] - M[1][2] * M[2][1];
    d12 = M[1][0] * M[2][2] - M[1][2] * M[2][0];
    d13 = M[1][0] * M[2][1] - M[1][1] * M[2][0];

    d21 = M[0][1] * M[2][2] - M[0][2] * M[2][1];
    d22 = M[0][0] * M[2][2] - M[0][2] * M[2][0];
    d23 = M[0][0] * M[2][1] - M[0][1] * M[2][0];

    d31 = M[0][1] * M[1][2] - M[0][2] * M[1][1];
    d32 = M[0][0] * M[1][2] - M[0][2] * M[1][0];
    d33 = M[0][0] * M[1][1] - M[0][1] * M[1][0];

    deti = 1.0 / (M[0][0] * d11 - M[0][1] * d12 + M[0][2] * d13);

    Inv(0, 0) =  deti * d11;
    Inv(0, 1) = -deti * d21;
    Inv(0, 2) =  deti * d31;
    Inv(1, 0) = -deti * d12;
    Inv(1, 1) =  deti * d22;
    Inv(1, 2) = -deti * d32;
    Inv(2, 0) =  deti * d13;
    Inv(2, 1) = -deti * d23;
    Inv(2, 2) =  deti * d33;

    return Inv;
}

inline
Vector3
operator * (const Matrix33& A, const Vector3& v)
{
    Vector3 Av = v;
    A.Multiply(Av);
    return Av;
}

inline
Matrix33
operator * (const Matrix33& A, const Matrix33& B)
{
    Matrix33 C;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            double c = 0.0;
            for (int l = 0; l < 3; ++l)
            {
                c += A(i, l) * B(l, j);
            }
            C(i, j) = c;
        }
    }
    return C;
}

inline
std::ostream&
operator << (std::ostream& out, const Matrix33& m)
{
    out << "[" << m(0, 0) << ", " << m(0, 1) << ", " << m(0, 2) << "]" << std::endl
        << "[" << m(1, 0) << ", " << m(1, 1) << ", " << m(1, 2) << "]" << std::endl
        << "[" << m(2, 0) << ", " << m(2, 1) << ", " << m(2, 2) << "]" << std::endl;
    return out;
}

#endif // INCLUDED_MATRIX_33_H__

