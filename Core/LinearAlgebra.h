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
// $Id: LinearAlgebra.h 116 2011-08-11 04:17:07Z kato $
#ifndef INCLUDED_LINEAR_ALGEBRA_H__
#define INCLUDED_LINEAR_ALGEBRA_H__

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <iostream>

namespace LinearAlgebra
{

template <class T, class E>
class VectorExpr
{
public:
    E Expr;
    VectorExpr(const E& expr) : Expr(expr) {}

    T At(int i) const { return Expr.At(i); }
    T operator [] (int i) const { return At(i); }
    int Size() const { return Expr.Size(); }
};

template <class T>
class Vector
{
public:
    Vector()
    : N(0), Data(NULL)
    {
    }

    Vector(const Vector<T>& v)
    : N(v.N), Data(new T[N])
    {
        SetTo(v);
    }

    Vector(int n)
    : N(n), Data(new T[N])
    {
    }

    ~Vector()
    {
        delete[] Data;
    }

    Vector<T>& Resize(int n)
    {
        delete[] Data;
        N = n;
        Data = new T[N];
        return *this;
    }

    Vector<T>& SetTo(const Vector<T>& v)
    {
        assert(N == v.N);
        for (int i = 0; i < N; ++i)
            Data[i] = v[i];
        return *this;
    }

    Vector<T>& SetTo(const T& v)
    {
        for (int i = 0; i < N; ++i)
            Data[i] = v;
        return *this;
    }

    int Size() const { return N; }

    const T& At(int i) const { return Data[i]; }
    T& At(int i) { return Data[i]; }
    const T& operator [] (int i) const { return Data[i]; }
    T& operator [] (int i) { return Data[i]; }
    const T& operator () (int i) const { return Data[i]; }
    T& operator () (int i) { return Data[i]; }

    //Vector<T>& operator = (T v) { SetTo(v); return *this; }
    Vector<T>& operator = (const T& v) { SetTo(v); return *this; }
    Vector<T>& operator = (const Vector<T>& v) { SetTo(v); return *this; }
    template <class E> Vector<T>& operator = (const VectorExpr<T, E>& e)
    {
        for (int i = 0; i < N; ++i)
            Data[i] = e.At(i);
        return *this;
    }

    Vector<T>& operator *= (const T& a)
    {
        for (int i = 0; i < N; ++i)
            Data[i] *= a;
        return *this;
    }

    T operator / (const T& a) const
    {
        Vector<T> res(N);
        for (int i = 0; i < N; ++i)
            res[i] = Data[i] / a;
        return res;
    }

    Vector<T>& operator /= (const T& a)
    {
        for (int i = 0; i < N; ++i)
            Data[i] /= a;
        return *this;
    }

protected:

private:
    int N;
    T* Data;
};

template <class T, class V>
class ScalarAddition
{
public:
    ScalarAddition(const T& a, const V& v) : mA(a), mV(v) {}
    T At(int i) const { return mA + mV.At(i); }
    int Size() const { return mV.Size(); }

private:
    T mA;
    const V& mV;
};

template <class T, class V1, class V2>
class VectorAddition
{
public:
    VectorAddition(const V1& v1, const V2& v2, T coef = T(1.0)) : mV1(v1), mV2(v2), mCoef(coef) {}
    T At(int i) const { return mV1.At(i) + mCoef * mV2.At(i); }
    int Size() const { return mV1.Size(); }

private:
    const V1& mV1;
    const V2& mV2;
    T mCoef;
};

template <class T, class V>
class ScalarMultiply
{
public:
    ScalarMultiply(const T& a, const V& v) : mA(a), mV(v) {}
    T At(int i) const { return mA * mV.At(i); }
    int Size() const { return mV.Size(); }

private:
    T mA;
    const V& mV;
};

template <class T, class V>
class Negate
{
public:
    Negate(const V& v) : mV(v) {}
    T At(int i) const { return -mV.At(i); }
    int Size() const { return mV.Size(); }

private:
    const V& mV;
};

template <class T>
VectorExpr<T, ScalarAddition<T, Vector<T> > >
operator + (const T& a, const Vector<T>& v)
{
    return VectorExpr<T, ScalarAddition<T, Vector<T> > >(ScalarAddition<T, Vector<T> >(a, v));
}

template <class T, class E>
VectorExpr<T, ScalarAddition<T, VectorExpr<T, E> > >
operator + (const T& a, const VectorExpr<T, E>& v)
{
    return VectorExpr<T, ScalarAddition<T, VectorExpr<T, E> > >(ScalarAddition<T, VectorExpr<T, E> >(a, v));
}

template <class T>
VectorExpr<T, VectorAddition<T, Vector<T>, Vector<T> > >
operator + (const Vector<T>& v1, const Vector<T>& v2)
{
    return VectorExpr<T, VectorAddition<T, Vector<T>, Vector<T> > >(VectorAddition<T, Vector<T>, Vector<T> >(v1, v2));
}

template <class T, class E>
VectorExpr<T, VectorAddition<T, Vector<T>, VectorExpr<T, E> > >
operator + (const Vector<T>& v1, const VectorExpr<T, E>& v2)
{
    return VectorExpr<T, VectorAddition<T, Vector<T>, VectorExpr<T, E> > >(
        VectorAddition<T, Vector<T>, VectorExpr<T, E> >(v1, v2)
        );
}

template <class T>
Vector<T>&
operator += (Vector<T>& v, const T& a)
{
    v = a + v;
    return v;
}

template <class T>
Vector<T>&
operator += (Vector<T>& v, const Vector<T>& a)
{
    v = v + a;
    return v;
}

template <class T, class E>
Vector<T>&
operator += (Vector<T>& v, const VectorExpr<T, E>& a)
{
    v = v + a;
    return v;
}

template <class T>
VectorExpr<T, ScalarAddition<T, Vector<T> > >
operator - (const Vector<T>& v, const T& a)
{
    return VectorExpr<T, ScalarAddition<T, Vector<T> > >(ScalarAddition<T, Vector<T> >(-a, v));
}

template <class T, class E>
VectorExpr<T, ScalarAddition<T, VectorExpr<T, E> > >
operator - (const VectorExpr<T, E>& v, const T& a)
{
    return VectorExpr<T, ScalarAddition<T, VectorExpr<T, E> > >(ScalarAddition<T, VectorExpr<T, E> >(-a, v));
}

template <class T>
VectorExpr<T, VectorAddition<T, Vector<T>, Vector<T> > >
operator - (const Vector<T>& v1, const Vector<T>& v2)
{
    return VectorExpr<T, VectorAddition<T, Vector<T>, Vector<T> > >(VectorAddition<T, Vector<T>, Vector<T> >(v1, v2, T(-1.0)));
}

template <class T, class E>
VectorExpr<T, VectorAddition<T, Vector<T>, VectorExpr<T, E> > >
operator - (const Vector<T>& v1, const VectorExpr<T, E>& v2)
{
    return VectorExpr<T, VectorAddition<T, Vector<T>, VectorExpr<T, E> > >(
        VectorAddition<T, Vector<T>, VectorExpr<T, E> >(v1, v2, T(-1.0))
        );
}

template <class T>
Vector<T>&
operator -= (Vector<T>& v, const T& a)
{
    v = v - a;
    return v;
}

template <class T>
Vector<T>&
operator -= (Vector<T>& v, const Vector<T>& a)
{
    v = v - a;
    return v;
}

template <class T, class E>
Vector<T>&
operator -= (Vector<T>& v, const VectorExpr<T, E>& a)
{
    v = v - a;
    return v;
}

template <class T>
VectorExpr<T, ScalarMultiply<T, Vector<T> > >
operator * (const T& a, const Vector<T>& v)
{
    return VectorExpr<T, ScalarMultiply<T, Vector<T> > >(ScalarMultiply<T, Vector<T> >(a, v));
}

template <class T, class E>
VectorExpr<T, ScalarMultiply<T, VectorExpr<T, E> > >
operator * (const T& a, const VectorExpr<T, E>& v)
{
    return VectorExpr<T, ScalarMultiply<T, VectorExpr<T, E> > >(ScalarMultiply<T, VectorExpr<T, E> >(a, v));
}

template <class T>
VectorExpr<T, ScalarMultiply<T, Vector<T> > >
operator * (const Vector<T>& v, const T& a)
{
    return VectorExpr<T, ScalarMultiply<T, Vector<T> > >(ScalarMultiply<T, Vector<T> >(a, v));
}

template <class T, class E>
VectorExpr<T, ScalarMultiply<T, VectorExpr<T, E> > >
operator * (const VectorExpr<T, E>& v, const T& a)
{
    return VectorExpr<T, ScalarMultiply<T, VectorExpr<T, E> > >(ScalarMultiply<T, VectorExpr<T, E> >(a, v));
}

template <class T>
T dot(const Vector<T>& v1, const Vector<T>& v2)
{
    assert(v1.Size() == v2.Size());

    T res = v1[0] * v2[0];
    for (int i = 1; i < v1.Size(); ++i)
        res += v1[i] * v2[i];
    return res;
}

template <class T>
T norm(const Vector<T>& v)
{
    T res = dot(v, v);
    return std::sqrt(res);
}

template <class T>
std::ostream& operator << (std::ostream& o, const Vector<T>& v)
{
    for (int i = 0; i < v.Size(); ++i)
    {
        o << ' ' << v[i];
    }
    return o;
}

template <class T, class E>
std::ostream& operator << (std::ostream& o, const VectorExpr<T, E>& v)
{
    for (int i = 0; i < v.Size(); ++i)
    {
        o << ' ' << v[i];
    }
    return o;
}

//------------------------------------------------------------------------------

template <class T>
class Matrix
{
public:
    Matrix()
    : M(0), N(0), Data(NULL)
    {
    }

    Matrix(const Matrix<T>& Mat)
    : M(Mat.M), N(Mat.N), Data(new T[M * N])
    {
        SetTo(Mat);
    }

    Matrix(int m, int n)
    : M(m), N(n), Data(new T[M * N])
    {
    }

    ~Matrix()
    {
        delete[] Data;
    }

    int Rows() const { return M; }
    int Columns() const { return N; }

    Matrix<T>& SetTo(const Matrix<T>& Mat)
    {
        assert(M == Mat.M && N == Mat.N);
        for (int j = 0; j < M; ++j)
            for (int i = 0; i < N; ++i)
                At(i, j) = Mat(i, j);
        return *this;
    }

    Matrix<T>& SetTo(const T& v)
    {
        for (int j = 0; j < M; ++j)
            for (int i = 0; i < N; ++i)
                At(i, j) = v;
        return *this;
    }

    const T& At(int i, int j) const { return Data[i * M + j]; }
    T& At(int i, int j) { return Data[i * M + j]; }
    const T& operator() (int i, int j) const { return At(i, j); }
    T& operator() (int i, int j) { return At(i, j); }

    Matrix<T>& operator = (const T& v) { SetTo(v); return *this; }
    Matrix<T>& operator = (const Matrix<T>& Mat) { SetTo(Mat); return *this; }

protected:

private:
    int M, N;
    T* Data;
};

template <class T>
Vector<T> operator * (const Matrix<T>& A, const Vector<T>& x)
{
    Vector<T> b(A.Rows());
    b = 0.0;
    for (int i = 0; i < A.Rows(); ++i)
        for (int j = 0; j < A.Columns(); ++j)
            b[i] += A(i, j) * x[j];
    return b;
}

template <class T>
std::ostream& operator << (std::ostream& o, const Matrix<T>& m)
{
    for (int i = 0; i < m.Rows(); ++i)
    {
        for (int j = 0; j < m.Columns(); ++j)
        {
            o << ' ' << m(i, j);
        }
        o << std::endl;
    }
    return o;
}

}

#endif // INCLUDED_LINEAR_ALGEBRA_H__
