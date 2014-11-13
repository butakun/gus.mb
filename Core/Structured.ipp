// $Id: Structured.ipp 312 2013-12-13 08:50:20Z kato $

#include <cmath>

template <class T>
Structured<T>&
Structured<T>::Add(const T& coeff, const Structured<T>& A)
{
    assert(Range.Shape() == A.Range.Shape() && NDOF == A.NDOF);

    T* src = A.Data;
    T* dest = Data;
    for (int k = Range.Start.K; k <= Range.End.K; ++k)
    {
        for (int j = Range.Start.J; j <= Range.End.J; ++j)
        {
            for (int i = Range.Start.I; i <= Range.End.I; ++i)
            {
                T* dst = At(i, j, k);
                T* src = A(i, j, k);
                for (int l = 0; l < NDOF; ++l)
                {
                    dst[l] += coeff * src[l];
                }
            }
        }
    }
    return *this;
}

template <class T>
Structured<T>&
Structured<T>::Add(const T& coeff, const Structured<T>& A, const IndexRange& range)
{
    assert(NDOF == A.NDOF);

    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                T* dst = At(i, j, k);
                T* src = A(i, j, k);
                for (int l = 0; l < NDOF; ++l)
                {
                    dst[l] += coeff * src[l];
                }
            }
        }
    }

    return *this;
}

template <class T>
Structured<T>&
Structured<T>::Multiply(const T& coeff)
{
    T* dest = Data;
    for (int k = 0; k < KDim; ++k)
    {
        for (int j = 0; j < JDim; ++j)
        {
            for (int i = 0; i < IDim; ++i)
            {
                for (int l = 0; l < NDOF; ++l)
                {
                    *dest = *dest * coeff;
                    ++dest;
                }
            }
        }
    }
    return *this;
}

template <class T>
Structured<T>&
Structured<T>::SetTo(const Structured<T>& A)
{
    assert(Range.Shape() == A.Range.Shape() && NDOF == A.NDOF);

    T* src = A.Data;
    T* dest = Data;
    for (int k = 0; k < KDim; ++k)
    {
        for (int j = 0; j < JDim; ++j)
        {
            for (int i = 0; i < IDim; ++i)
            {
                for (int l = 0; l < NDOF; ++l)
                {
                    *dest++ = *src++;
                }
            }
        }
    }
    return *this;
}

template <class T>
Structured<T>&
Structured<T>::SetTo(const T& value)
{
    T* dest = Data;
    for (int i = 0; i < IDim; ++i)
    {
        for (int j = 0; j < JDim; ++j)
        {
            for (int k = 0; k < KDim; ++k)
            {
                for (int l = 0; l < NDOF; ++l)
                {
                    *dest++ = value;
                }
            }
        }
    }
    return *this;
}

template <class T>
Structured<T>&
Structured<T>::SetTo(const T* values)
{
    // FIXME: check if the shape of A matches mine.
    T* dest = Data;
    for (int i = 0; i < IDim; ++i)
    {
        for (int j = 0; j < JDim; ++j)
        {
            for (int k = 0; k < KDim; ++k)
            {
                const T* value = values;
                for (int l = 0; l < NDOF; ++l)
                {
                    *dest++ = *value++;
                }
            }
        }
    }
    return *this;
}

template <class T>
Structured<T>&
Structured<T>::SetTo(const T* values, const IndexRange& range)
{
    for (int i = range.Start.I; i <= range.End.I; ++i)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int k = range.Start.K; k <= range.End.K; ++k)
            {
                const T* value = values;
                T* dest = At(i, j, k);
                for (int l = 0; l < NDOF; ++l)
                {
                    *dest = *value++;
                }
            }
        }
    }
    return *this;
}

template <class T>
T
Structured<T>::L2Norm(const IndexRange& range) const
{
    T res = 0.0;
    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                T* d = At(i, j, k);
                for (int l = 0; l < NDOF; ++l)
                {
                    res += d[l] * d[l];
                }
            }
        }
    }
    res = std::sqrt(res);
    return res;
}

template <class T>
void
Structured<T>::ReduceSquared(T* sq, T* sqMax, std::vector<IndexIJK>& maxlocs, const IndexRange& range) const
{
    for (int l = 0; l < NDOF; ++l)
    {
        sq[l] = 0.0;
    }
    for (int l = 0; l < NDOF + 1; ++l)
    {
        sqMax[l] = 0.0;
    }
    maxlocs.clear();
    maxlocs.resize(NDOF + 1);

    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                T* d = At(i, j, k);
                T sqAll = 0.0;
                for (int l = 0; l < NDOF; ++l)
                {
                    T s = d[l] * d[l];
                    sq[l] += s;
                    sqAll += s;
                    if (s > sqMax[l])
                    {
                        sqMax[l] = s;
                        maxlocs[l] = IndexIJK(i, j, k);
                    }
                }
                if (sqAll > sqMax[NDOF])
                {
                    sqMax[NDOF] = sqAll;
                    maxlocs[NDOF] = IndexIJK(i, j, k);
                }
            }
        }
    }
}

template <class T>
T
Structured<T>::RMS(const IndexRange& range) const
{
    T res = 0.0;
    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                T* d = At(i, j, k);
                for (int l = 0; l < NDOF; ++l)
                {
                    res += d[l] * d[l];
                }
            }
        }
    }
    res = std::sqrt(res / T(range.Count()));
    return res;
}

template <class T>
void
Structured<T>::RMS(T* rms, const IndexRange& range) const
{
    for (int l = 0; l < NDOF; ++l)
    {
        rms[l] = 0.0;
    }

    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                T* d = At(i, j, k);
                for (int l = 0; l < NDOF; ++l)
                {
                    rms[l] += d[l] * d[l];
                }
            }
        }
    }

    for (int l = 0; l < NDOF; ++l)
    {
        rms[l] = std::sqrt(rms[l] / T(range.Count()));
    }
}

template <class T>
void
Structured<T>::RMS(T* rms, IndexIJK* maxlocs, const IndexRange& range) const
{
    T* maxvals = new T[NDOF];
    for (int l = 0; l < NDOF; ++l)
    {
        rms[l] = 0.0;
        maxlocs[l] = IndexIJK(0, 0, 0);
        maxvals[l] = At(0, 0, 0)[l] * At(0, 0, 0)[l];
    }

    for (int k = range.Start.K; k <= range.End.K; ++k)
    {
        for (int j = range.Start.J; j <= range.End.J; ++j)
        {
            for (int i = range.Start.I; i <= range.End.I; ++i)
            {
                T* d = At(i, j, k);
                for (int l = 0; l < NDOF; ++l)
                {
                    T dd = d[l] * d[l];
                    rms[l] += dd;
                    if (dd > maxvals[l])
                    {
                        maxvals[l] = dd;
                        maxlocs[l] = IndexIJK(i, j, k);
                    }
                }
            }
        }
    }

    for (int l = 0; l < NDOF; ++l)
    {
        rms[l] = std::sqrt(rms[l] / T(range.Count()));
    }

    delete[] maxvals;
}

template <class T>
std::ostream&
Structured<T>::Dump(std::ostream& o, const IndexRange& r, const char* key) const
{
    for (int k = r.Start.K; k <= r.End.K; ++k)
    {
        for (int j = r.Start.J; j <= r.End.J; ++j)
        {
            for (int i = r.Start.I; i <= r.End.I; ++i)
            {
                o << IndexIJK(i, j, k) << " " << key << ":";
                for (int l = 0; l < NDOF; ++l)
                {
                    o << " " << At(i, j, k)[l];
                }
                o << std::endl;
            }
        }
    }
    return o;
}

template <class T>
Structured<T>&
Structured<T>::operator *= (const T& value)
{
    T* dest = Data;
    for (int i = 0; i < IDim; ++i)
    {
        for (int j = 0; j < JDim; ++j)
        {
            for (int k = 0; k < KDim; ++k)
            {
                for (int l = 0; l < NDOF; ++l)
                {
                    *dest = *dest * value;
                    ++dest;
                }
            }
        }
    }
    return *this;
}

#if 0
template<class T>
template<class E>
Structured<T>&
Structured<T>::operator = (const StructuredExpr<T, E>& expr)
{
    // FIXME
    expr.EvaluateAndAssign(*this);
    return *this;
}

template<class T, class E>
StructuredExpr<T, E>::StructuredExpr(const E& expr)
:   mExpr(expr)
{
}

template<class T, class E>
void
StructuredExpr<T, E>::EvaluateAndAssign(Structured<T>& S) const
{
    IndexRange r = S.GetRange();
    int dof = S.DOF();
    for (int k = r.Start.K; k <= r.End.K; ++k)
        for (int j = r.Start.J; j <= r.End.J; ++j)
            for (int i = r.Start.I; i <= r.End.I; ++i)
                for (int l = 0; l < dof; ++l)
                    S(i, j, k)[l] = mExpr.EvaluateAt(i, j, k, l);
}

#else
template<class T>
template<class E>
Structured<T>&
Structured<T>::operator = (const StructuredExpr<T, E>& expr)
{
    // FIXME
    //expr.EvaluateAndAssign(*this);

    IndexRange r = GetRange();
    int dof = DOF();
    for (int k = r.Start.K; k <= r.End.K; ++k)
        for (int j = r.Start.J; j <= r.End.J; ++j)
            for (int i = r.Start.I; i <= r.End.I; ++i)
                for (int l = 0; l < dof; ++l)
                    At(i, j, k, l) = expr.At(i, j, k, l);
    return *this;
}

template<class T>
template<class E>
Structured<T>&
Structured<T>::operator = (const E& expr)
{
    // FIXME
    //expr.EvaluateAndAssign(*this);

    IndexRange r = GetRange();
    int dof = DOF();
    for (int k = r.Start.K; k <= r.End.K; ++k)
        for (int j = r.Start.J; j <= r.End.J; ++j)
            for (int i = r.Start.I; i <= r.End.I; ++i)
                for (int l = 0; l < dof; ++l)
                    At(i, j, k, l) = expr.At(i, j, k, l);
    return *this;
}

template <class T>
template <class F>
void
Structured<T>::Apply(F& f)
{
    IndexRange r = GetRange();
    for (int k = r.Start.K; k <= r.End.K; ++k)
        for (int j = r.Start.J; j <= r.End.J; ++j)
            for (int i = r.Start.I; i <= r.End.I; ++i)
                f(At(IndexIJK(i, j, k)));
}

#endif
