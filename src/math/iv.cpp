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

#include "iv.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <boost/math/special_functions/next.hpp>

#include "helpers.h"

Iv::Iv(const double val)
    : interval(val), exactly_zero(false)
{
}

Iv::Iv( const double lower, const double upper)
    : interval(lower,upper), exactly_zero(false)
{
    check();
}

Iv::Iv( const bool not_exactly_zero)
    : interval(0), exactly_zero(true)
{
    if (not_exactly_zero)
        throw("Iv::Iv(true) not allowed, meant to construct exactly zero Iv!");
}

Iv::~Iv()
{
}

/*!
 * \param rhs Interval to add to this interval
 *
 * Calls operation in IvSSE2.
 * Also provides Iv operator+(const Iv&, const Iv&) due to boost::arithmetic.
 */
Iv &Iv::operator+= (const Iv &rhs)
{
    if (exactly_zero) {
        *this = rhs;
    } else if (~rhs.exactly_zero) {
        interval += rhs.interval;
        check();
    }
    return *this;
}

/*!
 * \param rhs Interval to subtract from this interval
 *
 * Calls operation in IvSSE2.
 * Also provides Iv operator-(const Iv&, const Iv&) due to boost::arithmetic.
 */
Iv &Iv::operator-= (const Iv &rhs)
{
    if (exactly_zero) {
        *this = -rhs;
    } else if (~rhs.exactly_zero) {
        interval -= rhs.interval;
        check();
    }
    return *this;
}

/*!
 * \param rhs Interval to multiply this interval by
 *
 * Calls operation in IvSSE2.
 * Also provides Iv operator*(const Iv&, const Iv&) due to boost::arithmetic.
 *
 * Note that multiplying an interval of exactly zero by another interval results in exactly zero, even if the
 * other interval is infinity! Our implementation assumes that a "zero interval" actually stands for a non-existing
 * term, e.g., the imaginary part of a complex interval in case this represents a real number!
 */
Iv &Iv::operator*= (const Iv &rhs)
{
    if (exactly_zero || rhs.exactly_zero) {
        *this = Iv(false);
    } else {
        interval *= rhs.interval;
        check();
    }

    return *this;
}

/*!
 * \param rhs Interval to divide this interval by
 *
 * Calls operation in IvSSE2.
 * Also provides Iv operator*(const Iv&, const Iv&) due to boost::arithmetic.
 *
 * Note that dividing an interval of exactly zero by another interval results in exactly zero (except the other interval
 * is zero as well), even if the other interval is infinity! Our implementation assumes that a "zero interval" actually
 * stands for a non-existing term, e.g., the imaginary part of a complex interval in case this represents a real number!
 */
Iv &Iv::operator/= (const Iv &rhs)
{
    if (exactly_zero && ~rhs.exactly_zero) {
        *this = Iv(false);
    } else {
        *this *= rhs.inv();
    }

    return *this;
}

/*!
 * Calls operation in IvSSE2.
 */
Iv Iv::operator- () const
{
    Iv ret;
    ret.interval = -interval;
    ret.exactly_zero = exactly_zero;

    return ret;
}

Iv Iv::inv() const
{
    Iv ret;
    ret.interval = 1.0/interval;

    return ret;
}

/*!
 * \param i Non-negative integer representing the exponent
 *
 * The power-function actually calculates the resulting interval using double precision arithmetics and then rounds the
 * lower bound downwards and the upper bound uppwards (except the lower bound if it is zero as the original interval
 * included zero)
 */
Iv Iv::power(unsigned int i) const
{
    Iv ret(1.0);
    for (size_t j=0; j<i; j++)
        ret *= *this;
    return ret;

    if (i==0) {
        return Iv(1.0);
    } else if ( (i%2==0) && (lower() < 0.0) ){ //even exponent with negative part
        if ( upper() < 0.0) { //even exponent, completely negative
            return Iv(rndD(pow(-upper(),(int)i)),rndU(pow(-lower(),(int)i)));
        } else { //even exponent, positive and negative
            if (upper() > -lower())
                return Iv(0.0,rndU(pow(upper(),(int)i)));
            else
                return Iv(0.0,rndU(pow(-lower(),(int)i)));
        }
    } else { // odd exponent with negative part or positive
        return Iv(rndD(pow(lower(),(int)i)),rndU(pow(upper(),(int)i)));
    }
}

/*!
 * Makes use of sincos().
 * \bug Calculation of m1 and m2 may not be safe, i.e., we may wrongly determine in which quadrant of the
 * sine/cosine we are.
 */
Iv Iv::sin() const
{
    if (width() >= 2*M_PI) // two extrema included
        return Iv(-1.0,1.0);
    else {
        double s1=std::sin(lower());
        double s2=std::sin(upper());

        double m1=std::fmod(lower()-M_PI_2,2*M_PI);
        if (m1<0.0)
            m1 += 2*M_PI;
        double m2=std::fmod(upper()-M_PI_2,2*M_PI);
        if (m2<0.0)
            m2 += 2*M_PI;

        //std::cout << "s1=" << s1 << ", s2=" << s2 << ", m1=" << m1 << ", m2=" << m2 << std::endl;
        return sincos(m1,m2,s1,s2);
    }
}

/*!
 * Makes use of sincos().
 * \bug Calculation of m1 and m2 may not be safe, i.e., we may wrongly determine in which quadrant of the
 * sine/cosine we are.
 */
Iv Iv::cos() const
{
    if (width() >= 2*M_PI) // two extrema included
        return Iv(-1.0,1.0);
    else {
        double s1=std::cos(lower());
        double s2=std::cos(upper());

        double m1=std::fmod(lower(),2*M_PI);
        if (m1<0.0)
            m1 += 2*M_PI;
        double m2=std::fmod(upper(),2*M_PI);
        if (m2<0.0)
            m2 += 2*M_PI;

        //std::cout << "s1=" << s1 << ", s2=" << s2 << ", m1=" << m1 << ", m2=" << m2 << std::endl;
        return sincos(m1,m2,s1,s2);
    }
}

Iv Iv::exp() const
{
    double a,b;

    if (lower() == 0.0)
        a=1.0;
    else if (lower() == -INFINITY)
        a = 0.0;
    else
        a = rndD(std::exp(lower()));

    if (upper() == 0.0)
        b = 1.0;
    else
        b = rndU(std::exp(upper()));

    return Iv(a,b);
}

/*!
 * \param m1 Quadrant of lower bound (m1>M_PI => rising, m1<=M_PI => falling)
 * \param m2 Quadrant of upper bound (m1>M_PI => rising, m1<=M_PI => falling)
 * \param s1 Sine/cosine of lower bound (double precision, not rounded outwards yet)
 * \param s2 Sine/cosine of upper bound (double precision, not rounded outwards yet)
 * \return Interval of sine/cosine
 *
 * This function checks whether there is an extremum of the sine/cosine between lower and upper bound and determines the
 * overall result of the trigonometric function accordingly.
 *
 * The non-exact values (not 1.0 or -1.0) are rounded outwards to account for rounding errors floating-point arithmetics.
 *
 * \warning If m1 and m2 are not necessarily correct this can lead to wrong results!
 */
Iv Iv::sincos(double m1, double m2, double s1, double s2) const
{
    if (width() > M_PI) {    // one or two extrema included
        if (m1 > M_PI) {     // lower bound on rising part
            if (m2 > M_PI)   // upper bound on rising part
                return Iv(-1.0,1.0);
            else             // upper bound on falling part
                return Iv(rndD(s1<s2?s1:s2),1.0);
        } else {             // lower bound on falling part
            if (m2 > M_PI)   // upper bound on rising part
                return Iv(-1.0,rndU(s1>s2?s1:s2));
            else             // upper bound on falling part
                return Iv(-1.0,1.0);
        }
    } else {                 // no or one extremum
        if (m1 > M_PI) {     // lower bound on rising part
            if (m2 > M_PI)   // upper bound on rising part
                return Iv(rndD(s1),rndU(s2));
            else             // upper bound on falling part
                return Iv(rndD(s1<s2?s1:s2),1.0);
        } else {             // lower bound on falling part
            if (m2 > M_PI)   // upper bound on rising part
                return Iv(-1.0,rndU(s1>s2?s1:s2));
            else             // upper bound on falling part
                return Iv(rndD(s2),rndU(s1));
        }
    }
}

/*!
 * \param factor Double value to scale th interval by
 * \warning This function is not safe, i.e., no outward rounding is applied.
 */
Iv Iv::scale(double factor) const
{
    double rad = width() / 2.0 * factor;
    return Iv(mid()-rad,mid()+rad);
}

Iv &Iv::operator= (const Iv &rhs)
{
    interval = rhs.interval;
    exactly_zero = rhs.exactly_zero;

    return *this;
}

Vec<Iv> Iv::splitAt(const Vec<double> &split_points) const
{
    if (split_points.empty()) {
        fprintf(stderr, "Error in Iv::splitAt: split_points empty!\n");
        exit(1);
    }

    Vec<double> limits = lower();

    limits.push_back(split_points);
    limits.push_back(upper());

    Vec<Iv> parts;
    for (size_t i = 0; i < limits.size()-1; i++) {
        if (limits[i+1] < limits[i]) {
            fprintf(stderr, "Error in Iv::splitAt: limits %d smaler than limits %d!\n",(int)(i+1),(int)i);
            exit(1);
        }
        parts.push_back(Iv(limits[i],limits[i+1]));
    }

    return parts;
}

Vec<Iv> Iv::splitAt(double split_point) const
{
    Vec<double> split_points(split_point);
    return splitAt(split_points);
}

Vec<Iv> Iv::splitEqual(size_t num_parts) const
{
    Vec<double> split_points;

    for (size_t i = 1; i < num_parts; i++)
        split_points.push_back(lower() + i*(width()/num_parts));

    return splitAt(split_points);
}

Iv Iv::intersect(const Iv &i) const
{
    if (isNaN() || i.isNaN()) {
        //printf("WARINING: INTERSECT 1 or 2 is NaN!\n");
        return Iv(NAN);
    }

    double b = lower();
    double t = upper();

    if (i.lower() > b)
        b = i.lower();

    if (i.upper() < t)
        t = i.upper();

    return Iv(b,t);
}

Iv Iv::unify(const Iv &i) const
{
    if (isNaN() || i.isNaN()) {
        //printf("WARINING: UNION 1 or 2 is NaN!\n");
        return Iv(NAN);
    }

    double b = lower();
    double t = upper();

    if (i.lower() < b)
        b = i.lower();

    if (i.upper() > t)
        t = i.upper();

    return Iv(b,t);
}

bool Iv::includes(double val) const
{
    return ( (lower() <= val) && (upper() >= val));
}

bool Iv::includes(const Iv &val) const
{
    return ( (lower() <= val.lower()) && (upper() >= val.upper()) );
}

double Iv::boundViolation(double x) const
{
    double du = upper() - x;
    double dl = x - lower();

    if ( (du>=0.0) && (dl>=0.0) )
        return (du>dl?dl:du);
    else if ( (du>=0.0) || (dl>=0.0) )
        return (du<dl?du:dl);
    else
        return (du<dl?dl:du);
}

void Iv::print() const
{
    std::cout << interval << std::endl;
}

bool Iv::isNaN() const
{
    return (isnan(lower()) || isnan(upper()));
}

bool Iv::isFinite() const
{
    return (std::isfinite(lower()) && std::isfinite(upper()));
}

void Iv::check()
{
    if (lower() > upper()) {
        //std::cerr << "Iv::check(): lower() > upper()!" << std::endl;
        interval = IvSSE2<>(NAN,NAN);
    }
}

double Iv::rndU(double x) const
{
    if (!std::isfinite(x))
        return x;
    else
        return boost::math::float_next(x);
}

double Iv::rndD(double x) const
{
    if (!std::isfinite(x))
        return x;
    else
        return boost::math::float_prior(x);
}

Iv sin(const Iv &iv)
{
    return iv.sin();
}

Iv cos(const Iv &iv)
{
    return iv.cos();
}

Iv exp(const Iv &iv)
{
    return iv.exp();
}

std::ostream & operator<<(std::ostream &os, const Iv& iv)
{
    os << iv.interval;

    return os;
}
