// $Id: Structured.h 310 2013-12-10 07:59:40Z kato $
#ifndef INCLUDED_STRUCTURED_H__
#define INCLUDED_STRUCTURED_H__

#include <vector>
#include <iostream>
#include <cassert>

class IndexIJK
{
public:
    int I, J, K;

    IndexIJK() : I(0), J(0), K(0) {}
    IndexIJK(int i, int j, int k) : I(i), J(j), K(k) {}
    IndexIJK(int i) : I(i), J(i), K(i) {}

    IndexIJK operator + (const IndexIJK& i2) const { return IndexIJK(I + i2.I, J + i2.J, K + i2.K); }
    IndexIJK operator - (const IndexIJK& i2) const { return IndexIJK(I - i2.I, J - i2.J, K - i2.K); }
    IndexIJK& operator += (const IndexIJK& i) { I += i.I; J += i.J; K += i.K; return *this; }
    IndexIJK& operator -= (const IndexIJK& i) { I -= i.I; J -= i.J; K -= i.K; return *this; }

    IndexIJK operator * (int a) const { return IndexIJK(I * a, J * a, K * a); }

    int operator [] (int i) const { return i == 0 ? I : (i == 1 ? J : K); }
    int& operator [] (int i) { return i == 0 ? I : (i == 1 ? J : K); }

protected:

private:
};

class IndexRange
{
public:
    // Start and End indices are inclusive, i.e. (0, 0) - (1, 0) is meant to cover (0, 0) and (1, 0).
    IndexIJK Start;
    IndexIJK End;

    IndexRange() : Start(0, 0, 0), End(0, 0, 0) {}
    IndexRange(const IndexRange& range) : Start(range.Start), End(range.End) {}
    IndexRange(const IndexIJK& start, const IndexIJK& end) : Start(start), End(end) {}
    IndexRange(int is, int js, int ks, int ie, int je, int ke) : Start(is, js, ks), End(ie, je, ke) {}

    int Count() const { IndexIJK d = Size(); return d.I * d.J * d.K; }
    IndexIJK Shape() const { return End + IndexIJK(1, 1, 1) - Start; }
    IndexIJK Size() const { return End + IndexIJK(1, 1, 1) - Start; }

    IndexRange Canonical() const;
    bool IsCanonical() const { return Start.I <= End.I && Start.J <= End.J && Start.K <= End.K; }
    IndexRange& Canonize() { *this = Canonical(); return *this; }

    IndexRange& operator = (const IndexRange& range) { Start = range.Start; End = range.End; return *this; }

protected:

private:
};

inline
std::ostream& operator << (std::ostream& o, const IndexIJK& ijk)
{
    o << '(' << ijk.I << ',' << ijk.J << ',' << ijk.K << ')';
    return o;
}

inline
std::ostream& operator << (std::ostream& o, const IndexRange& range)
{
    o << range.Start << " - " << range.End;
    return o;
}

inline
bool operator == (const IndexIJK& i1, const IndexIJK& i2)
{
    return i1.I == i2.I && i1.J == i2.J && i1.K == i2.K;
}

inline
bool operator == (const IndexRange& r1, const IndexRange& r2)
{
    return r1.Start == r2.Start && r1.End == r2.End;
}

inline
IndexRange operator + (const IndexRange& r, const IndexIJK& i)
{
    return IndexRange(r.Start + i, r.End + i);
}

class IndexIterator
{
public:
    IndexIterator(const IndexRange& range)
    : mRange(range), mIndex(mRange.Start)
    {
        assert(mRange.IsCanonical());
    }
    const IndexIJK& Index() const { return mIndex; }
    void Advance()
    {
        ++mIndex.I;
        if (mIndex.I > mRange.End.I)
        {
            mIndex.I = mRange.Start.I;
            ++mIndex.J;
        }
        if (mIndex.J > mRange.End.J)
        {
            mIndex.J = mRange.Start.J;
            ++mIndex.K;
        }
    }
    bool IsEnd() const
    { 
        return mIndex.K > mRange.End.K;
    }

protected:

private:
    IndexRange mRange;
    IndexIJK mIndex;
};

template<class T> class Structured;

template<class T, class E>
class StructuredExpr
{
public:
    StructuredExpr(const E& expr) : mExpr(expr) {}
    T At(int i, int j, int k, int l) const { return mExpr.At(i, j, k, l); }

protected:

private:
    E mExpr;
};

template <class T>
class Structured
{
public:
    T* Data;

    Structured()
    :   Data(NULL), NDOF(0), Range(0, 0, 0, 0, 0, 0), IDim(0), JDim(0), KDim(0)
    {
    }

    Structured(const Structured<T>& s)
    :   Data(s.Data), NDOF(s.NDOF), Range(s.Range), IDim(s.IDim), JDim(s.JDim), KDim(s.KDim)
    {
    }

    Structured(T* data, int ndof, const IndexRange& range)
    :   Data(data), NDOF(ndof), Range(range)
    {
        IndexIJK size = range.Size();
        IDim = size.I;
        JDim = size.J;
        KDim = size.K;
    }

    // Allocates the data, (but the caller is responsible for freeing it)
    Structured(int ndof, const IndexRange& range)
    :   Data(new T[range.Count() * ndof]), NDOF(ndof), Range(range)
    {
        IndexIJK size = range.Size();
        IDim = size.I;
        JDim = size.J;
        KDim = size.K;
    }

    void Allocate(int ndof, const IndexRange& range)
    {
        delete[] Data;
        Data = new T[range.Count() * ndof];
        NDOF = ndof;
        Range = range;
        IndexIJK size = range.Size();
        IDim = size.I;
        JDim = size.J;
        KDim = size.K;
    }

    const IndexRange& GetRange() const { return Range; }
    void SetRange(const IndexRange& range)
    {
        Range = range;
        IndexIJK size = range.Size();
        IDim = size.I; JDim = size.J; KDim = size.K;
    }

    //T* At(const IndexIJK& ijk) const { IndexIJK ii = ijk - Range.Start; return &Data[(ii.I * JDim * KDim + ii.J * KDim + ii.K) * NDOF]; }
    T* At(const IndexIJK& ijk) const
    {
        IndexIJK ii = ijk - Range.Start;
        return &Data[(ii.K * JDim * IDim + ii.J * IDim + ii.I) * NDOF];
    }
    T* At(int i, int j, int k) const { return At(IndexIJK(i, j, k)); }
    T* operator () (const IndexIJK& ijk) const { return At(ijk); }
    T* operator () (int i, int j, int k) const { return At(i, j, k); }
    const T& At(int i, int j, int k, int l) const { return At(i, j, k)[l]; }
    T& At(int i, int j, int k, int l) { return At(i, j, k)[l]; }

    Structured<T>& SetTo(const T& value);
    Structured<T>& SetTo(const Structured<T>& A);
    Structured<T>& SetTo(const T* values); // values is an array of NDOF elements.
    Structured<T>& SetTo(const T* values, const IndexRange& range); // values is an array of NDOF elements.
    Structured<T>& Add(const T& coeff, const Structured<T>& A);
    Structured<T>& Add(const T& coeff, const Structured<T>& A, const IndexRange& range);
    Structured<T>& Multiply(const T& coeff);

    template <class F> void Apply(F& f);

    int DOF() const {return this->NDOF; }
    void SetDOF(int dof) { this->NDOF = dof; }

    T L2Norm(const IndexRange& range) const;
    void ReduceSquared(T* sq, T* sqMax, std::vector<IndexIJK>& maxlocs, const IndexRange& range) const;
    T RMS(const IndexRange& range) const;
    void RMS(T* rms, const IndexRange& range) const;
    void RMS(T* rms, IndexIJK* maxlocs, const IndexRange& range) const;

    Structured<T>& operator = (const Structured<T>& A) { SetTo(A); return *this; }
    Structured<T>& operator = (const T& value) { SetTo(value); return *this; }
    Structured<T>& operator *= (const T& value);

#if 0
    template<class E> Structured<T>& operator = (const StructuredExpr<T, E>& expr);
#else
    template<class E> Structured<T>& operator = (const StructuredExpr<T, E>& expr);
    template<class E> Structured<T>& operator = (const E& expr);
#endif

    std::ostream& Dump(std::ostream& o, const IndexRange& r, const char* key) const;

protected:

private:
    int NDOF;
    IndexRange Range;
    int IDim, JDim, KDim;
};

template<class T>
std::ostream&
operator << (std::ostream& o, const Structured<T>& S)
{
    S.Dump(o, S.GetRange(), "");
    return o;
}

template<class T, class S>
class ScalarAddition
{
public:
    ScalarAddition(const T& a, const S& s)
    :   mA(a), mS(s)
    {
    }

    T At(int i, int j, int k, int l) const
    {
        return mA + mS.At(i, j, k, l);
    }

private:
    T mA;
    const S& mS;
};

template<class T, class S>
class ScalarMultiply
{
public:
    ScalarMultiply(const T& a, const S& s)
    :   mA(a), mS(s)
    {
    }

    T At(int i, int j, int k, int l) const
    {
        return mA * mS.At(i, j, k, l);
    }

private:
    T mA;
    const S& mS;
};

template<class T>
StructuredExpr<T, ScalarAddition<T, Structured<T> > >
operator + (const T& a, const Structured<T>& s)
{
    return StructuredExpr<T, ScalarAddition<T, Structured<T> > >(ScalarAddition<T, Structured<T> >(a, s));
}

template<class T, class E>
StructuredExpr<T, ScalarAddition<T, StructuredExpr<T, E> > >
operator + (const T& a, const StructuredExpr<T, E>& s)
{
    return StructuredExpr<T, ScalarAddition<T, StructuredExpr<T, E> > >(ScalarAddition<T, StructuredExpr<T, E> >(a, s));
}

template<class T>
StructuredExpr<T, ScalarAddition<T, Structured<T> > >
operator + (const Structured<T>& s, const T& a)
{
    return StructuredExpr<T, ScalarAddition<T, Structured<T> > >(ScalarAddition<T, Structured<T> >(a, s));
}

template<class T, class E>
StructuredExpr<T, ScalarAddition<T, StructuredExpr<T, E> > >
operator + (const StructuredExpr<T, E>& s, const T& a)
{
    return StructuredExpr<T, ScalarAddition<T, StructuredExpr<T, E> > >(ScalarAddition<T, StructuredExpr<T, E> >(a, s));
}

template<class T>
StructuredExpr<T, ScalarMultiply<T, Structured<T> > >
operator * (const T& a, const Structured<T>& s)
{
    return StructuredExpr<T, ScalarMultiply<T, Structured<T> > >(ScalarMultiply<T, Structured<T> >(a, s));
}

template<class T, class E>
StructuredExpr<T, ScalarMultiply<T, StructuredExpr<T, E> > >
operator * (const T& a, const StructuredExpr<T, E>& s)
{
    return StructuredExpr<T, ScalarMultiply<T, StructuredExpr<T, E> > >(ScalarMultiply<T, StructuredExpr<T, E> >(a, s));
}

template<class T>
StructuredExpr<T, ScalarMultiply<T, Structured<T> > >
operator * (const Structured<T>& s, const T& a)
{
    return StructuredExpr<T, ScalarMultiply<T, Structured<T> > >(ScalarMultiply<T, Structured<T> >(a, s));
}

template<class T, class E>
StructuredExpr<T, ScalarMultiply<T, StructuredExpr<T, E> > >
operator * (const StructuredExpr<T, E>& s, const T& a)
{
    return StructuredExpr<T, ScalarMultiply<T, StructuredExpr<T, E> > >(ScalarMultiply<T, StructuredExpr<T, E> >(a, s));
}

template<class T>
StructuredExpr<T, ScalarMultiply<T, Structured<T> > >
operator / (const Structured<T>& s, const T& a)
{
    return StructuredExpr<T, ScalarMultiply<T, Structured<T> > >(ScalarMultiply<T, Structured<T> >(T(1.0) / a, s));
}

template<class T, class E>
StructuredExpr<T, ScalarMultiply<T, StructuredExpr<T, E> > >
operator / (const StructuredExpr<T, E>& s, const T& a)
{
    return StructuredExpr<T, ScalarMultiply<T, StructuredExpr<T, E> > >(ScalarMultiply<T, StructuredExpr<T, E> >(T(1.0) / a, s));
}

#include "Structured.ipp"

#endif // INCLUDED_STRUCTURED_H__

