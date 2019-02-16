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

#ifndef IV_SSE2_ROUNDING_H
#define IV_SSE2_ROUNDING_H

#include <emmintrin.h>

/*!
 * \brief Switch rounding mode for interval arithmetics
 *
 * Constructing this object with template parameter true (or empty) switches rounding to upward and stores old setting
 * of rounding, deconstruction then resets rounding.
 *
 * Constructing and destroying this object with template argument false does nothing!
 */
template <bool switch_rounding=true>
class IvSSE2rounding
{
public:
    IvSSE2rounding()
    {
        /* Save rounding modes */
        msr_orig = _mm_getcsr();

        /* Set to upward rounding */
        _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    }

    ~IvSSE2rounding()
    {
        /* Restore them */
        _mm_setcsr(msr_orig);
    }

private:
    unsigned msr_orig;
};

/*!
 * \brief Do not switch rounding mode for interval arithmetics
 *
 * Constructing and destroying this object with template argument false does nothing!
 */
template<>
class IvSSE2rounding<false>
{
public:
    IvSSE2rounding() {}
    ~IvSSE2rounding() {}
};

#endif // IV_SSE2_ROUNDING_H
