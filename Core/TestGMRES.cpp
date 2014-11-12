// $Id: TestGMRES.cpp 116 2011-08-11 04:17:07Z kato $

#include "LinearAlgebra.h"
#include "GMRES.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace LinearAlgebra;

#if 0
class BlockVector : public Vector<double>
{
public:
    BlockVector() : Vector<double>() {}
    BlockVector(int n) : Vector<double>(n) {}
    BlockVector(const BlockVector& v) : Vector<double>(v) {}
    BlockVector(const Vector<double>& v) : Vector<double>(v) {}
    ~BlockVector() {}

    BlockVector& AXPBY(double a, double b, const BlockVector& y)
    {
        for (int i = 0; i < Size(); ++i)
        {
            At(i) = a * At(i) + b * y[i];
        }
        return *this;
    }

    BlockVector& operator = (double v) { SetTo(v); return *this; }

protected:

private:
};
#endif

class NoPreconditioner
{
public:
    NoPreconditioner() {}
    void Solve(Vector<double>& res, const Vector<double>& b) const
    {
        for (int i = 0; i < b.Size(); ++i)
            res[i] = b[i];
    }
};

class JacobiPreconditioner
{
public:
    JacobiPreconditioner(const Matrix<double>& A_, int loop = 20)
    :   A(A_), x(A.Columns()), mLoop(loop)
    {
    }

    void Solve(Vector<double>& res, const Vector<double>& b)
    {
        res = 0.0;

        for (int loop = 0; loop < mLoop; ++loop)
        {
            for (int i = 0; i < b.Size(); ++i)
            {
                double axk = 0.0;
                for (int j = 0; j < b.Size(); ++j)
                {
                    if (j == i) continue;
                    axk += A(i, j) * res[j];
                }

                x[i] = (-axk + b[i]) / A(i, i);
            }
            res = x;
        }
    }

private:
    const Matrix<double>& A;
    Vector<double> x;
    int mLoop;
};

#if 0
class BlockVector : public Vector<double>
{
public:
    BlockVector(int n) : Vector<double>(n) {}

protected:

private:
};
#endif

class Operator
{
public:
    Operator(const Matrix<double>& A_) : A(A_) {}

    Vector<double>& InitializeBigVector(Vector<double>& v) const { v.Resize(A.Rows()); return v; }

    Vector<double> operator * (const Vector<double>& x) const { return A * x; }

    void Apply(Vector<double>& res, const Vector<double>& x) const { res = A * x; }

private:
    const Matrix<double>& A;
};

void GenerateMatrix(Matrix<double>& A)
{
    int m = A.Rows(), n = A.Columns();

    for (int i = 0; i < m; ++i)
    {
#if 0
        for (int j = 0; j < n; ++j)
        {
            A(i, j) = 10.0 * double(std::rand()) / double(RAND_MAX);
        }
#else
        double d, e, f;
        d = double(std::rand()) / double(RAND_MAX);
        e = double(std::rand()) / double(RAND_MAX);
        f = double(std::rand()) / double(RAND_MAX);
        d = d + e + f;

        A(i, i) = 2.0 * d;
        if (i > 0) A(i, i - 1) = e;
        if (i < (m - 1)) A(i, i + 1) = f;
#endif
    }
}

int
main()
{
    // Random seed
    std::srand(int(time(NULL)));

    const int N = 15;
    Vector<double> v1(N);

    v1 = 1.234;
    std::cout << v1 << std::endl;

    v1 += 2.0;
    std::cout << v1 << std::endl;

    std::cout << norm(v1) << std::endl;

    Matrix<double> m1(N, N);
    m1 = 1.234;
    std::cout << m1 << std::endl;

    Matrix<double> A(N, N);
    GenerateMatrix(A);
    std::cout << "A = " << std::endl << A << std::endl;

    Vector<double> b(N);
    for (int i = 0; i < N; ++i)
        b[i] = 10.0 * double(std::rand()) / double(RAND_MAX);
    std::cout << "b = " << std::endl << b << std::endl;

    JacobiPreconditioner M(A, 1);
    NoPreconditioner M2;
    Operator Ax(A);

    //GMRES<Operator, Vector<double>, JacobiPreconditioner, Vector<double>, Matrix<double>, double> gmres(12, Ax, M);
    GMRES<Operator, Vector<double>, NoPreconditioner, Vector<double>, Matrix<double>, double> gmres(12, Ax, M2);

    Vector<double> x(N), x2(N);
    x = 0.0;
    //x.SetTo(0.0);

#if 0
    M.Solve(x2, b);
    std::cout << "b = " << std::endl << b << std::endl;
    std::cout << "Ax = " << std::endl << A * x2 << std::endl;
    std::cout << "Ax - b = " << std::endl << A * x2 - b << std::endl;
    return 0;
#endif

    int iter;
    double tol = 0.01;
    gmres.Solve(x, b, iter, tol);

    std::cout << "iter = " << iter << std::endl;

    std::cout << "X = " << std::endl << x << std::endl;

    Vector<double> Axmb(N);
    Axmb = A * x;
    std::cout << "b = " << std::endl << b << std::endl;
    std::cout << "Ax = " << std::endl << Axmb << std::endl;
    Axmb = A * x - b;
    std::cout << "Ax - b = " << std::endl << Axmb << std::endl;

    return 0;
}

