// $Id: Reconstructor.h 254 2012-10-07 17:51:59Z kato $
#ifndef INCLUDED_RECONSTRUCTOR_H__
#define INCLUDED_RECONSTRUCTOR_H__

template <int N, class LIMITER, class CODEC>
class Reconstructor
{
public:
    Reconstructor();

    void Reconstruct(
        double* UL, double* UR,
        const Structured<double>& U,
        const Structured<double>& Rad,
        const IndexIJK& I, const IndexIJK& dI
        ) const;

protected:

private:
    LIMITER mLimiter;
    CODEC mCodec;
};

template <int N>
class FirstOrderReconstructor
{
public:
    FirstOrderReconstructor();

    void Reconstruct(
        double* UL, double* UR,
        const Structured<double>& U,
        const Structured<double>& Rad,
        const IndexIJK& I, const IndexIJK& dI
        ) const;

protected:

private:
};

#if 0
inline
double MinModLimiter(double r)
{
    return std::max(0.0, std::min(1.0, r));
}
#else
class MinModLimiter
{
public:
    MinModLimiter() {}
    double operator () (double r) const
    {
        return std::max(0.0, std::min(1.0, r));
    }
};
#endif

class PrimitiveVariableCodec
{
public:
    PrimitiveVariableCodec();

    // Conservative (U) -> Primitive (Q)
    void Decode(double* Q, double* U, double* Rad) const;

    // Primitive (Q) -> Conservative (U)
    void Encode(double* U, double* Q, double* Rad) const;

protected:

private:
    double mGamma;
};

#include "Reconstructor.ipp"

#endif // INCLUDED_RECONSTRUCTOR_H__

