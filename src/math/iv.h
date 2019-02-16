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

#ifndef IV_H
#define IV_H

#include "vec.hpp"
#include "iv_sse2.h"
#include <cmath>
#include <boost/operators.hpp>

/*!
 * \brief The Iv class implements intervals.
 *
 * This includes arithmetic operations, some intrinsic functions (only the ones needed for this project) as well as
 * several other functions on intervals.
 *
 * The actual implementation of the basic arithmetic operations is done in IvSSE2 and requires a processor with SSE2
 * insturctionsis.
 *
 * \todo It would be interesting to replace this implementation of interval arithmetics with some existing library
 * (e.g, boost interval or CGAL). However, most available widely spread implementations (e.g., the ones named) do not
 * contain any optimization for SSE2 as used here. I am not sure how big the overall performance impact of this would
 * be.
 */
class Iv : boost::arithmetic<Iv>
{
public:
    Iv(const double val = NAN);                 //!< Construct interval with lower and upper bound set to val.
    Iv(const double lower, const double upper); //!< Construct interval with given lower and upper bound.
    Iv(const bool not_exactly_zero);            //!< Construct interval which is exactly zero using Iv(false)
    virtual ~Iv();                              //!< Destructor.

    Iv &operator= (const Iv &rhs);   //!< Assignment operator
    Iv &operator+= (const Iv &rhs);  //!< Addition.
    Iv &operator-= (const Iv &rhs);  //!< Subtraction.
    Iv &operator*= (const Iv &rhs);  //!< Multiplication.
    Iv &operator/= (const Iv &rhs);  //!< Division.
    Iv operator-() const;            //!< Negation

    Iv inv() const;                  //!< Inverse 1.0/Iv
    Iv power(unsigned int i) const;  //!< Integer power
    Iv sin() const;                  //!< Sine
    Iv cos() const;                  //!< Cosine
    Iv exp() const;                  //!< Exponential

    Iv scale(double factor=0) const; //!< Scale interval by a given factor (mid-point stays the same)

    bool isNaN() const;              //!< Is upper or lower bound NaN? \todo Should be called NaI?
    bool isFinite() const;           //!< Is upper and lower bound finite?

    void print() const;              //!< Print interval to std::cout

    double lower() const {return interval.lower();}            //!< Return lower bound
    double upper() const {return interval.upper();}            //!< Return upper bound
    double mid() const   {return ((lower()+upper()) * 0.5);}   //!< Return midpoint as double \warning Rounding errors possible
    double width() const {return (upper()-lower());}           //!< Return width as double \warning Rounding errors possible

    /*!
     * \brief Split interval into several intervals at given points
     * \param split_points Vector of n points at which the interval is split, must be sorted in asscending order and lie
     * within the interval
     * \returns Vector of n+1 intervals
     */
    Vec<Iv> splitAt(const Vec<double> &split_points) const;
    /*!
     * \brief Split interval into two intervals at given points
     * \param split_point Point at which the interval is split, must lie within the interval
     * \returns Vector of two intervals
     */
    Vec<Iv> splitAt(double split_point) const;
    /*!
     * \brief Split interval into a given number of intervals with equal width
     * \param num_parts Number of intervals to return
     * \returns Vector of num_parts intervals
     * \warning The split points are determined using floating point arithmetics, i.e., the resulting interals may
     * differe in width minimally.
     * \bug The upper bound of the new intervals must not necessarily be identical to the upper bound of the original
     * interval! Use this function with caution or fix!
     */
    Vec<Iv> splitEqual(size_t num_parts) const;

    //! Upper bound of interval less than lower bound of rhs
    bool operator<  (const Iv &rhs) const {return (upper() <  rhs.lower());}
    //! Upper bound of interval less than or equal to lower bound of rhs
    bool operator<= (const Iv &rhs) const {return (upper() <= rhs.lower());}
    //! Lower bound of interval greater than upper bound of rhs
    bool operator>  (const Iv &rhs) const {return (lower() >  rhs.upper());}
    //! Lower bound of interval greater than or equal to upper bound of rhs
    bool operator>= (const Iv &rhs) const {return (lower() >= rhs.upper());}
    //! Lower and upper bounds both identical
    bool operator== (const Iv &rhs) const {return ((lower() == rhs.lower())&&(upper() == rhs.upper()));}
    //! Lower and upper bounds not both identical
    bool operator!= (const Iv &rhs) const {return ((lower() != rhs.lower())||(upper() != rhs.upper()));}

    /*!
     * \brief Intersection of two intervals
     * \param i second interval
     * \return Intersection, NaN-Interval if the two intervals do not intersect.
     */
    Iv intersect(const Iv &i) const;

    /*!
     * \brief Union of two intervals
     * \param i second interval
     * \return Union as interval, if there is a gap between the two intervals this is contained in the result.
     */
    Iv unify(const Iv &i) const;

    /*!
     * \brief Checks whether the interval includes a double value
     * \param val Value to check for inclusion
     */
    bool includes(double val) const;
    /*!
     * \brief Checks whether the interval completely includes another interval
     * \param val Interval to check for inclusion
     */
    bool includes(const Iv &val) const;
    /*!
     * \brief Checks whether a value is contained in the interval and returns its location relative to the bounds
     * \param x Value for which inclusion is checked
     * \return Positive value indicates x is contained in interval with value indicating distance from closer boundary,
     * negative value indicates x is not contained in interval and -x is distance to closer boundary.
     */
    double boundViolation(double x) const;

protected:
    //! Check whether upper bound is greater or equal to lower bound, otherwise set to NaN-interval
    void check();

    //! Helper function used to calculate sine and cosine
    Iv sincos(double m1, double m2, double s1, double s2) const;

    double rndU(double x) const; //!< Return x rounded upwards (smallest representable number which is larger than x)
    double rndD(double x) const; //!< Return x rounded downwards (largest representable number which is smaller than x)

    IvSSE2<> interval; //!< Actual interval value
    /*!
     * \brief Interval is exactly zero.
     * \warning Exactly zero has the meaning of a non-existing term, e.g., the imaginary part of a complex interval in
     * case this represents a real number! This implies, e.g., that multiplying an interval by an interval which is
     * exactly zero results in an interval which is exactly zero, even if the other interval is infinite or NaN/NaI.
     */
    bool exactly_zero;

    //! Stream operator for debug-output
    friend std::ostream & operator<<(std::ostream &os, const Iv& iv);
};

//! Calculate sine of interval using Iv::sin()
Iv sin(const Iv &iv);
//! Calculate cosine of interval using Iv::cos()
Iv cos(const Iv &iv);
//! Calculate exponential of interval using Iv::exp()
Iv exp(const Iv &iv);

//! Stream operator for debug-output
std::ostream & operator<<(std::ostream &os, const Iv& iv);

#endif // IV_H
