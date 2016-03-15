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
#ifndef INCLUDE_GMRES_H__
#define INCLUDE_GMRES_H__
// $Id: GMRES.h 300 2013-09-28 09:04:50Z kato $
//
// GMRES solver template, based on the netlib template available from
// http://www.netlib.org/templates/
// modified by hiromasa.kato@cenaero.be
// further modified by hiromasa@iwate-u.ac.jp

//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the 
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************


template <class Matrix, class Vector, class BigVector>
void 
Update(BigVector &x, int k, Matrix &h, Vector &s, BigVector v[])
{
    Vector y(s);

    // Backsolve:  
    for (int i = k; i >= 0; i--)
    {
        y(i) /= h(i,i);
        for (int j = i - 1; j >= 0; j--)
        {
            y(j) -= h(j,i) * y(i);
        }
    }

    for (int j = 0; j <= k; j++)
    {
        x += v[j] * y(j);
    }
}

template < class Real >
Real 
abs(Real x)
{
    return (x > 0 ? x : -x);
}

#include <cmath> 

template<class Real> 
void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
    if (dy == 0.0)
    {
        cs = 1.0;
        sn = 0.0;
    }
    else if (abs(dy) > abs(dx))
    {
        Real temp = dx / dy;
        sn = 1.0 / sqrt( 1.0 + temp*temp );
        cs = temp * sn;
    }
    else
    {
        Real temp = dy / dx;
        cs = 1.0 / sqrt( 1.0 + temp*temp );
        sn = temp * cs;
    }
}

template<class Real> 
void ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
    Real temp  =  cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
}

template <class Operator, class BigVector, class Preconditioner, class Vector, class Matrix, class Real>
class GMRES
{
public:
    GMRES(
        int m_,
        Operator& A_,
        Preconditioner& M_
        );
    virtual ~GMRES();

    virtual int Solve(BigVector& x, const BigVector& b, int& max_iter, Real& tol);

protected:

private:
    int m;
    Operator& A;
    Preconditioner& M;

    // workspace
    BigVector r;
    BigVector* v;
};

template <class Operator, class BigVector, class Preconditioner, class Vector, class Matrix, class Real>
GMRES<Operator, BigVector, Preconditioner, Vector, Matrix, Real>::GMRES(
    int m_,
    Operator& A_,
    Preconditioner& M_
    )
:   m(m_), A(A_), M(M_)
{
    A.InitializeBigVector(r);
    v = new BigVector[m + 1];
    for (int i = 0; i < m + 1; ++i)
    {
        A.InitializeBigVector(v[i]);
    }
}

template <class Operator, class BigVector, class Preconditioner, class Vector, class Matrix, class Real>
GMRES<Operator, BigVector, Preconditioner, Vector, Matrix, Real>::
~GMRES()
{
    delete[] v;
    v = NULL;
}

template <class Operator, class BigVector, class Preconditioner, class Vector, class Matrix, class Real>
int 
GMRES<Operator, BigVector, Preconditioner, Vector, Matrix, Real>::
Solve(BigVector& x, const BigVector& b, int& max_iter, Real& tol)
{
    Real resid;
    int i, j = 1, k;
    Vector s(m + 1), cs(m + 1), sn(m + 1);
    Matrix H(m + 1, m);
    BigVector w, Ax;
    A.InitializeBigVector(w);
    A.InitializeBigVector(Ax);
  
    max_iter = m;

    //Real normb = norm(M.Solve(b));
    M.Solve(r, b); // M * r = b;
    Real normb = norm(r);

    //-----------------------------------------
    //r = M.Solve(b - A * x);
    //-----------------------------------------
    A.Apply(Ax, x); // Ax = A * x
    //Ax.AXPBY(-1.0, 1.0, b);
    Ax = b - Ax;
    M.Solve(r, Ax); // M * r = Ax
    Real beta = norm(r);
  
    if (normb == 0.0)
    {
        normb = 1;
    }
  
    //if ((resid = norm(r) / normb) <= tol)
    resid = beta / normb;
    std::cout << "resid = " << resid << std::endl;
    if (resid <= tol)
    {
        tol = resid;
        max_iter = 0;
        return 0;
    }

    while (j <= max_iter)
    {
        v[0] = r * (1.0 / beta);    // ??? r / beta
        s = 0.0;
        s(0) = beta;
    
        for (i = 0; i < m && j <= max_iter; i++, j++)
        {
            //w = M.Solve(A * v[i]);
            A.Apply(Ax, v[i]);
            M.Solve(w, Ax);
            for (k = 0; k <= i; k++)
            {
                H(k, i) = dot(w, v[k]);
                w -= H(k, i) * v[k];
            }
            H(i + 1, i) = norm(w);
            v[i + 1] = w * (1.0 / H(i + 1, i)); // ??? w / H(i+1, i)
  
            for (k = 0; k < i; k++)
            {
                ApplyPlaneRotation(H(k, i), H(k + 1, i), cs(k), sn(k));
            }
            
            GeneratePlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
            ApplyPlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
            ApplyPlaneRotation(s(i), s(i + 1), cs(i), sn(i));
            
            resid = abs(s(i + 1)) / normb;
            std::cout << j << ", resid = " << resid << std::endl;
            if (resid < tol)
            {
                Update(x, i, H, s, v);
                tol = resid;
                max_iter = j;
                return 0;
            }
        }
        Update(x, m - 1, H, s, v);
        //---------------------------------
        //r = M.Solve(b - A * x);
        A.Apply(Ax, x); // Ax = A * x
        //Ax.AXPBY(-1.0, 1.0, b);
        Ax = b - Ax;
        M.Solve(r, Ax); // M * r = Ax
        //---------------------------------
        beta = norm(r);
        resid = beta / normb;
        std::cout << j << ", resid = " << resid << std::endl;
        if (resid < tol)
        {
            tol = resid;
            max_iter = j;
            return 0;
        }
    }
  
    tol = resid;

    return 1;
}

#endif // INCLUDE_GMRES_H__

