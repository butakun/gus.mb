// $Id: Vector3.h 228 2012-05-07 10:46:59Z kato $
#ifndef INCLUDED_VECTOR3_H__
#define INCLUDED_VECTOR3_H__

#include <iostream>
#include <cmath>

class Vector3
{
public:
                Vector3() {}
                Vector3(double x, double y, double z) { mXYZ[0] = x; mXYZ[1] = y; mXYZ[2] = z; }
                Vector3(const Vector3& v) { mXYZ[0] = v.X(); mXYZ[1] = v.Y(); mXYZ[2] = v.Z(); }
                Vector3(const Vector3& v1, const Vector3& v2) { mXYZ[0] = v2.X() - v1.X(); mXYZ[1] = v2.Y() - v1.Y(); mXYZ[2] = v2.Z() - v1.Z(); }
                Vector3(double* v) { mXYZ[0] = v[0]; mXYZ[1] = v[1]; mXYZ[2] = v[2]; }
                Vector3(float* v) { mXYZ[0] = v[0]; mXYZ[1] = v[1]; mXYZ[2] = v[2]; }
                Vector3(double* v1, double* v2) { mXYZ[0] = v2[0] - v1[0]; mXYZ[1] = v2[1] - v1[1]; mXYZ[2] = v2[2] - v1[2]; }

    const double* Data() const { return mXYZ; }
    double      X() const { return mXYZ[0]; }
    double      Y() const { return mXYZ[1]; }
    double      Z() const { return mXYZ[2]; }
#ifndef SWIG
    double&     X() { return mXYZ[0]; }
    double&     Y() { return mXYZ[1]; }
    double&     Z() { return mXYZ[2]; }
#endif

    double      MagSq() const { return X() * X() + Y() * Y() + Z() * Z(); }
    double      Mag() const { return std::sqrt(MagSq()); }

    void        Normalize() { double a = 1.0 / Mag(); mXYZ[0] *= a; mXYZ[1] *= a; mXYZ[2] *= a; }
    Vector3     Normalized() const { Vector3 n(*this); n.Normalize(); return n; }

    void        GetValues(double* v) const { v[0] = mXYZ[0]; v[1] = mXYZ[1]; v[2] = mXYZ[2]; }

    void        Apply(double M[][3]);
    static void RotationMatrix(double R[3][3], const Vector3& rotation);

#ifndef SWIG
    bool        operator == (const Vector3& v) const { return mXYZ[0] == v[0] && mXYZ[1] == v[1] && mXYZ[2] == v[2]; }
    bool        operator != (const Vector3& v) const { return mXYZ[0] != v[0] || mXYZ[1] != v[1] || mXYZ[2] != v[2]; }
    double      operator [] (size_t i) const { return mXYZ[i]; }
    double&     operator [] (size_t i) { return mXYZ[i]; }
    void        operator = (const Vector3& v) { mXYZ[0] = v.X(); mXYZ[1] = v.Y(); mXYZ[2] = v.Z(); }
    Vector3     operator - () const { return Vector3(-X(), -Y(), -Z()); }
    Vector3&    operator += (const double a) { mXYZ[0] += a; mXYZ[1] += a; mXYZ[2] += a; return *this; }
    Vector3&    operator -= (const double a) { mXYZ[0] -= a; mXYZ[1] -= a; mXYZ[2] -= a; return *this; }
    Vector3&    operator *= (const double a) { mXYZ[0] *= a; mXYZ[1] *= a; mXYZ[2] *= a; return *this; }
    Vector3&    operator /= (const double a) { mXYZ[0] /= a; mXYZ[1] /= a; mXYZ[2] /= a; return *this; }
#endif

protected:

private:
    double      mXYZ[3];
};

inline
void
Vector3::Apply(double M[][3])
{
    double x2, y2, z2;
    x2 = M[0][0] * X() + M[0][1] * Y() + M[0][2] * Z();
    y2 = M[1][0] * X() + M[1][1] * Y() + M[1][2] * Z();
    z2 = M[2][0] * X() + M[2][1] * Y() + M[2][2] * Z();
    mXYZ[0] = x2;
    mXYZ[1] = y2;
    mXYZ[2] = z2;
}

inline
void
Vector3::RotationMatrix(double R[3][3], const Vector3& rot)
{
    Vector3 u(0.0, 0.0, 0.0);
    if (rot != Vector3(0.0, 0.0, 0.0))
    {
        u = rot.Normalized();
    }
    double theta = rot.Mag();
    double cos = std::cos(theta), sin = std::sin(theta), onemcos = 1.0 - cos;

    R[0][0] = cos + u.X() * u.X() * onemcos;
    R[1][0] = u.Y() * u.X() * onemcos + u.Z() * sin;
    R[2][0] = u.Z() * u.X() * onemcos - u.Y() * sin;

    R[0][1] = u.X() * u.Y() * onemcos - u.Z() * sin;
    R[1][1] = cos + u.Y() * u.Y() * onemcos;
    R[2][1] = u.Z() * u.Y() * onemcos + u.X() * sin;

    R[0][2] = u.X() * u.Z() * onemcos + u.Y() * sin;
    R[1][2] = u.Y() * u.Z() * onemcos - u.X() * sin;
    R[2][2] = cos + u.Z() * u.Z() * onemcos;
}

inline
double
dot_product(const Vector3& v1, const Vector3& v2)
{
    return v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z();
}

inline
Vector3
cross_product(const Vector3& v1, const Vector3& v2)
{
    return Vector3(
            v1.Y() * v2.Z() - v1.Z() * v2.Y(),
            v1.Z() * v2.X() - v1.X() * v2.Z(),
            v1.X() * v2.Y() - v1.Y() * v2.X()
            );
}

inline
Vector3
operator * (double a, const Vector3& v)
{
    return Vector3(a * v.X(), a * v.Y(), a * v.Z());
}

inline
Vector3
operator * (const Vector3& v, double a)
{
    return Vector3(a * v.X(), a * v.Y(), a * v.Z());
}

inline
Vector3
operator / (const Vector3& v, double a)
{
    return Vector3(v.X() / a, v.Y() / a, v.Z() / a);
}

inline
Vector3
operator + (const Vector3& v1, const Vector3& v2)
{
    return Vector3(v1.X() + v2.X(), v1.Y() + v2.Y(), v1.Z() + v2.Z());
}

inline
Vector3
operator - (const Vector3& v1, const Vector3& v2)
{
    return Vector3(v1.X() - v2.X(), v1.Y() - v2.Y(), v1.Z() - v2.Z());
}

inline
std::ostream&
operator << (std::ostream& out, const Vector3& v)
{
    out << "(" << v.X() << ", " << v.Y() << ", " << v.Z() << ")";
    return out;
}

#endif // INCLUDED_VECTOR3_H__

