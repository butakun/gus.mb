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
std::ostream&
operator << (std::ostream& out, const Matrix33& m)
{
    out << "[" << m(0, 0) << ", " << m(0, 1) << ", " << m(0, 2) << "]" << std::endl
        << "[" << m(1, 0) << ", " << m(1, 1) << ", " << m(1, 2) << "]" << std::endl
        << "[" << m(2, 0) << ", " << m(2, 1) << ", " << m(2, 2) << "]" << std::endl;
    return out;
}

#endif // INCLUDED_MATRIX_33_H__

