/*
    Copyright (C) 2017 Thomas Schauss

    This file is part of glob_stab.

    glob_stab is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    glob_stab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with glob_stab. If not, see <http://www.gnu.org/licenses/>.
*/

#include "iv_sse2.h"

template <bool switch_rounding>
IvSSE2<switch_rounding>::IvSSE2(__m128d iv)
{
    val = iv;
}

/*!
 * Store negative lower bound -> all computations using rounding to +infty
 */
template <bool switch_rounding>
IvSSE2<switch_rounding>::IvSSE2(double lowerupper)
{
    //
    val = (__m128d){-lowerupper, lowerupper};
}

/*!
 * Store negative lower bound -> all computations using rounding to +infty
 */
template <bool switch_rounding>
IvSSE2<switch_rounding>::IvSSE2(double lower, double upper)
{
    val = (__m128d){-lower, upper};
}

template <bool switch_rounding>
IvSSE2<switch_rounding>::~IvSSE2()
{

}

template <bool switch_rounding>
double IvSSE2<switch_rounding>::lower() const
{
    return -val_d[0];
}

template <bool switch_rounding>
double IvSSE2<switch_rounding>::upper() const
{
    return val_d[1];
}

template <bool switch_rounding>
IvSSE2<switch_rounding> &IvSSE2<switch_rounding>::operator=  (const IvSSE2 &rhs)
{
    val = rhs.val;
    return *this;
}

template <bool switch_rounding>
IvSSE2<switch_rounding> &IvSSE2<switch_rounding>::operator+= (const IvSSE2 &rhs)
{
    IvSSE2rounding<switch_rounding> rnd;
    val = add(val,rhs.val);
    return *this;
}

template <bool switch_rounding>
IvSSE2<switch_rounding> &IvSSE2<switch_rounding>::operator-= (const IvSSE2 &rhs)
{
    IvSSE2rounding<switch_rounding> rnd;
    val = add(val,neg(rhs.val));
    return *this;
}

template <bool switch_rounding>
IvSSE2<switch_rounding> &IvSSE2<switch_rounding>::operator*= (const IvSSE2 &rhs)
{
    IvSSE2rounding<switch_rounding> rnd;
    val = mul(val,rhs.val);
    return *this;
}

template <bool switch_rounding>
IvSSE2<switch_rounding> &IvSSE2<switch_rounding>::operator/= (const IvSSE2 &rhs)
{
    IvSSE2rounding<switch_rounding> rnd;

    /* make +0.0 to -0.0 -> correct result in inv(...) */
    IvSSE2 r(rhs);
    if (r.val_d[0] == 0.0)
        r.val_d[0] = -0.0;
    if (r.val_d[1] == 0.0)
        r.val_d[1] = -0.0;

    val = mul(val,inv(r.val));
    return *this;
}

template <bool switch_rounding>
IvSSE2<switch_rounding> IvSSE2<switch_rounding>::operator-() const
{
    return IvSSE2(neg(val));
}

/*!
 * Component-wise addition with rounding towards +infty results in correct interval addition.
 * Assumes rounding mode is towards positive infinity.
 * Taken from http://locklessinc.com/articles/interval_arithmetic/
 */
template <bool switch_rounding>
__m128d IvSSE2<switch_rounding>::add(__m128d x, __m128d y) const
{
    return x+y;
}

/*!
 * Swapping the high and low part of the 128-bit value results in an exact negation.
 * Taken from http://locklessinc.com/articles/interval_arithmetic/
 */
template <bool switch_rounding>
__m128d IvSSE2<switch_rounding>::neg(__m128d x) const
{
    return _mm_shuffle_pd(x, x, 1);
}

/*!
 * Multiplication is more complex bit can also be implemented efficiently using SSE2-commands.
 * Assumes rounding mode is towards positive infinity.
 * Taken from http://locklessinc.com/articles/interval_arithmetic/
 */
template <bool switch_rounding>
__m128d IvSSE2<switch_rounding>::mul(__m128d x, __m128d y) const
{
    __m128d t1 = (__m128d)_mm_shuffle_epi32((__m128i) x, 0xee);
    __m128d t2 = (__m128d)_mm_shuffle_epi32((__m128i) y, 0xee);

    __m128d t3 = _mm_xor_pd(x, t1);
    __m128d t4 = _mm_xor_pd(y, t2);

    if (_mm_movemask_pd(_mm_and_pd(t3, t4))) {
        __m128d c = {0.0, 0.0};
        __m128d c1 = _mm_cmple_pd(t2, c);
        __m128d c2 = _mm_cmple_pd(t1, c);
        __m128d c3 = (__m128d) {-0.0, 0.0};

        x = cSwap(_mm_xor_pd(x, c3), c1);
        y = cSwap(_mm_xor_pd(y, c3), c2);

        return x * _mm_xor_pd(y, c3);
    }

    /* There is a zero overlap */
    t1 = (__m128d)_mm_shuffle_epi32((__m128i) x, 0x4e) * _mm_unpacklo_pd(y, y);
    t2 *= x;

    return _mm_max_pd(t1, t2);
}

/*!
 * Calculation of the reciprocal is quite simple and efficient.
 * Assumes rounding mode is towards positive infinity.
 * Taken from http://locklessinc.com/articles/interval_arithmetic/
 */
template <bool switch_rounding>
__m128d IvSSE2<switch_rounding>::inv(__m128d x) const
{
    /* if x[0] < 0 and x[1] < 0  -> Iv contains 0 (changed from <=) */
    if (_mm_movemask_pd(_mm_cmplt_pd((__m128d) {0, 0}, x)) == 3) {
        return (__m128d){INFINITY, INFINITY};
    }
    x = _mm_shuffle_pd(x, x, 1);

    return _mm_div_pd((__m128d) {-1.0, -1.0}, x);
}

/*!
 * Taken from http://locklessinc.com/articles/interval_arithmetic/
 */
template <bool switch_rounding>
__m128d IvSSE2<switch_rounding>::cSwap(__m128d x, __m128d cond) const
{
    __m128d t = _mm_and_pd((__m128d)_mm_shuffle_epi32((__m128i) x, 0x4e), cond);

    cond = _mm_andnot_pd(cond, x);

    return _mm_or_pd(cond, t);
}

template class IvSSE2<false>;
template class IvSSE2<true>;

template <bool switch_rounding>
std::ostream & operator<<(std::ostream &os, const IvSSE2<switch_rounding>& iv)
{
    os << "[ " << iv.lower() << ", " << iv.upper() << " ]";

    return os;
}

template std::ostream & operator<<(std::ostream &os, const IvSSE2<true>& iv);
template std::ostream & operator<<(std::ostream &os, const IvSSE2<false>& iv);
