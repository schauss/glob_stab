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

#ifndef IV_SSE2_H
#define IV_SSE2_H

#include <cmath>
#include <emmintrin.h>
#include <iostream>
#include <boost/operators.hpp>

#include "iv_sse2_rounding.hpp"
/*!
 * \brief SSE2 Implementation of interval arithmetics
 *
 * Two modes are available:
 * * If consturcted using IvSSE<> or IvSSE<true> rounding mode is switched for arithmetic operations.
 * * If constructed using IvSSE<false> rounding mode is not switched. Then the user must take care of this (e.g. by
 *   constructiong IvSSE2rounding<> or IvSSE2rounding<true> which is also used to switch rounding mode here.
 *
 * The implementation of the actual arithmetics is based on http://locklessinc.com/articles/interval_arithmetic/ from
 * where most of the code was directly copied. The general idea is to store an interval in one 128-bit SSE2 register
 * with negative lower bound as this makes it possible to implement many operations very efficiently.
 *
 * This in turn is based on "INTERVAL ARITHMETIC USING SSE-2" by BRANIMIR LAMBOV which can be found at
 * http://www.daimi.au.dk/~barnie/intervals.pdf  and describes the implementation used in RealLib which is available
 * at https://github.com/blambov/RealLib and published under the Apache 2.0 license.
 *
 * \todo We could directly use RealLib for low-level interval calculations and keep the interface in Iv.
 */
template <bool switch_rounding=true>
class IvSSE2 : boost::arithmetic<IvSSE2<switch_rounding> >
{
public:
    IvSSE2(__m128d iv);                           //!< Copy constructor.
    IvSSE2(double lowerupper);                    //!< Set lower and upper bound to lowerupper.
    IvSSE2(double lower, double upper);           //!< Set lower bound to lower and upper bound to upper.
    ~IvSSE2();                                    //!< Destructor.

    double lower() const;                         //!< Return lower bound (-val_d[0])
    double upper() const;                         //!< Return upper bound (val_d[1])

    IvSSE2 &operator=  (const IvSSE2 &rhs);       //!< Assignment operator
    IvSSE2 &operator+= (const IvSSE2 &rhs);       //!< Switch rounding, add two intervals using add()
    IvSSE2 &operator-= (const IvSSE2 &rhs);       //!< Switch rounding, subtract two intervals using add() and neg()
    IvSSE2 &operator*= (const IvSSE2 &rhs);       //!< Switch rounding, multiply two intervals using mul()
    IvSSE2 &operator/= (const IvSSE2 &rhs);       //!< Switch rounding, divide two intervals using mul() and inv()
    IvSSE2 operator-() const;                     //!< Negate an interval using neg()

private:
    __m128d add(__m128d x, __m128d y) const;      //!< Sum of two intervals
    __m128d neg(__m128d x) const;                 //!< Negation of interval
    __m128d mul(__m128d x, __m128d y) const;      //!< Product of two intervals
    __m128d inv(__m128d x) const;                 //!< Inverse (or reciprocal) of an interval

    __m128d cSwap(__m128d x, __m128d cond) const; //!< Conditional swap using mask

    union {
        __m128d val;
        double val_d[2];                          //!< Negative lower bound (val_d[0]) and positive upper bound (val_d[1])
    };
};

//! Stream operator for debug-output
template <bool switch_rounding>
std::ostream & operator<<(std::ostream &os, const IvSSE2<switch_rounding>& iv);

#endif // IV_SSE2_H
