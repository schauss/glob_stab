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

#ifndef TAYLORMODEL_H
#define TAYLORMODEL_H

#include <iostream>
#include <string>
#include <boost/operators.hpp>
#include <complex>

#include "config.h"
#include "vec.hpp"
#include "polynomial.h"
#include "iv.h"

template <class T>
class TaylorModel;

//! Stream operator for debug-output
template <class T>
std::ostream & operator<<(std::ostream &os, const TaylorModel<T>& t);

/*!
 * \brief This class implements Taylor Models which consist of a multivariate
 * Polynomial and an interval remainder.
 *
 * Arithmetic operations, some intrinsic functions, and different bounding-
 * methods are implemented. This allows to efficiently evaluate arbitrary
 * functions (not all intrinsic functions are implemented).
 *
 * The maximum total order of the multivariate Polynomial (i.e., the maximum sum
 * of the values of the exponents of all Monomials) is limited to
 * max_total_degree. When evaluating intrinsic functions or arithmetic
 * operations, bounds are determined for Monomials above this degree and these
 * are added to the interval remainder.
 *
 * The method is described in \cite Berz1998. Details, e.g., on intrinsic
 * functions, can be found in \cite Makino2003.
 *
 * In contrast to those publications we do not automatically enlarge the
 * interval remainder to account for floating-point rounding-errors. Therefore,
 * TaylorModels of type double do not necessarily overapproximate the exact
 * solution as the interval remainder is purely determined from terms of an
 * order higher than the Taylor Model.
 *
 * To make sure the TaylorModel is an overapproximation of the exact value set
 * is you must use a TaylorModel of type interval (Iv) instead. In this case, a
 * polynomial with intervals as coefficients is used, and the foating-point
 * errors in the polynomial part are therefore determined for each term. Note
 * that this is compiutationally less efficient than the method proposed in
 * \cite Berz1998 and \cite Makino2003.
 *
 * \todo We should remove direct calls to boost::program_options and instead
 * make all configuration-values accessible using member functions!
 * \todo Many more intrinsic functions could be implemented.
 */
template <class T=double>
class TaylorModel
        : boost::arithmetic<TaylorModel<T>
        , boost::arithmetic2<TaylorModel<T>, double
        > >
{
public:
    //! Configuration options (see init and the class Config for more information)
    static boost::program_options::options_description options();
    /*!
     * \brief Initialize Taylor Model
     * \param var_names Specify variable names (sets dimension of Polynomial)
     *
     * The following configuration options are set when calling init:
     * * tm_order: maximum total degree of Taylor Model (default: 10). Higher
     * values lead to longer computation times but may approximate actual
     * function more exactly.
     * * tm_center: use domain [-1,1] instead of [0,1] ([0,1] when not set).
     * A centered Taylor expansion results in a smaller interval remainder but
     * requires a coordinate transformation if Berstein polynomials are used
     * to evaluate boundsimproves the ")
     * * tm_b_off: turn off bernstein bounding. If this is set, then naiive
     * interval arithmetics are used to determine bounds of the Taylor Model
     * (when bounds are requested explicitly and when bounding higher-order
     * terms). Otherwise the convex hull of the Bernstein coefficients is used.
     *
     * This function also calls Monomial::fixVariables(var_names) and thereby
     * fixes the variables available for multivariate Polynomials (and thereby
     * Taylor Models).
     */
    static void init(const Vec<std::string> &var_names);

    //! Return value of correct_trig. See correct_trig for more details.
    static bool correctTrigonometric() {return correct_trig;}
    //! Set value of correct_trig. See correct_trig for more details.
    static void correctTrigonometric(bool val) {correct_trig = val;}

    //! Construct TaylorModel which is exactly zero (polynomial: empty, remainder: empty)
    TaylorModel(bool not_exactly_zero = false);
    //! Construct TaylorModel consisting of one variable named var (remainder: empty)
    TaylorModel(const std::string& var);
    //! Construct TaylorModel consisting of one variable named var (remainder: empty)
    TaylorModel(const char var[]);
    //! Construct TaylorModel consisting of double value val (remainder: empty)
    TaylorModel(double val);
    //! Construct TaylorModel consisting of interval remainder (polynomial: empty)
    TaylorModel(const Iv &iv);
    //! Construct TaylorModel representing variable var on domain i (transform to configured domain)
    TaylorModel(const std::string &var, const Iv &i);
    //! Destructor.
    ~TaylorModel();

    //! Return Polynomial part
    const Polynomial<T> &polynomial() const {return poly;}
    /*!
     * \brief Return Polynomial part transformed to the domain [0,1] if
     * necessary.
     *
     * If a transformation is necessary the result is cached in poly_01.
     */
    const Polynomial<T> &polynomial_01() const;

    /*!
     * \brief Outer approximation of TaylorModel
     * \return Sum of outer approximation of polynomial and interval remainder.
     *
     * Depending on the setting of doBernsteinBound() in Polynomial this is
     * either evaluated using a Bernstein transformation or naiive interval
     * arithmetics.
     */
    Iv bound() const;
    //! Return remainder interval.
    const Iv &restBound() const {return rest;}

    /*!
     * \brief Return direction of maximum partial derivative.
     *
     * \todo Not implemented!
     */
    size_t maxDxDir() const {throw ("TaylorModel::maxDxDir not implemented!"); return 0;}

    //! Add rhs TaylorModel to this TaylorModel.
    TaylorModel<T> &operator+= (const TaylorModel<T> &rhs);
    //! Subtract rhs TaylorModel to this TaylorModel.
    TaylorModel<T> &operator-= (const TaylorModel<T> &rhs);
    //! Multiply this TaylorModel with rhs TaylorModel.
    TaylorModel<T> &operator*= (const TaylorModel<T> &rhs);
    //! Divide this TaylorModel by rhs TaylorModel.
    TaylorModel<T> &operator/= (const TaylorModel<T> &rhs);

    //! Add rhs double to this TaylorModel.
    TaylorModel<T> &operator+= (double rhs);
    //! Subtract rhs double to this TaylorModel.
    TaylorModel<T> &operator-= (double rhs);
    //! Multiply this TaylorModel with rhs double.
    TaylorModel<T> &operator*= (double rhs);
    //! Divide this TaylorModel by rhs double.
    TaylorModel<T> &operator/= (double rhs);

    //! Negation of TaylorModel
    TaylorModel<T> operator-() const;
    //! Reciprocal of TaylorModel
    TaylorModel<T> inv() const;
    //! Sine of TaylorModel
    TaylorModel<T> sin() const;
    //! Cosine of TaylorModel
    TaylorModel<T> cos() const;
    //! Exponential of TaylorModel.
    TaylorModel<T> exp() const;
    //! Return TaylorModel to the power i.
    TaylorModel<T> power(size_t i) const;

    //! Check whether the lower and upper bound of the TaylorModel both are identical to the double rhs-
    bool operator==(double rhs) const;

private:
    //! Return constant part of polynomial.
    Iv cF() const;
    //! Horner scheme used to evaluate sin, cos, exp, and inv.
    TaylorModel<T> horner(const TaylorModel<T> &x, const Vec<Iv> &a) const;
    //! Limit maximum total order of Polynomial to max_total_degree and add interval bound of higher orders to interval remainder.
    void cutoff();
    //! Polynomial part of TaylorModel
    Polynomial<T> poly;
    //! Cache for polynomial_01() (only used if domain != [0,1]).
    mutable Polynomial<T> poly_01;
    //! Remainder interval
    Iv rest;

    //! Domain of all variables (generally either [0,1] or [-1,1]
    static Iv domain;
    //! Maximum total order of Polynomial (i.e., maximum sum of exponents per Monomial)
    static int max_total_degree;
    /*!
     * \brief Replace TaylorModel of sin()/cos() with sin()/cos() of interval-
     * bound in some cases.
     *
     * If sin()/cos() of the Taylor Model results in NaN or the remainder
     * interval is larger than sin()/cos() of the interval-bound this
     * replacement is carried out. Note that sin()/cos() of the interval is
     * intersected with the interval [-1, 1] as this is the range of values
     * these functions can take.
     */
    static bool correct_trig;

    //! Stream operator for debug-output
    friend std::ostream & operator<< <>(std::ostream &os, const TaylorModel& t);
};

//! Free sin() function to allow writting sin(TaylorModel). Calls TaylorModel::sin()
template <class T> TaylorModel<T> sin(const TaylorModel<T> &t) {return t.sin();}
//! Free cos() function to allow writting cos(TaylorModel). Calls TaylorModel::cos()
template <class T> TaylorModel<T> cos(const TaylorModel<T> &t) {return t.cos();}
//! Free exp() function to allow writting exp(TaylorModel). Calls TaylorModel::exp()
template <class T> TaylorModel<T> exp(const TaylorModel<T> &t) {return t.exp();}
//! Exponential function of complex Taylor Model.
template <class T> std::complex<TaylorModel<T> > exp(const std::complex<TaylorModel<T> > &t);

//! Free sinh() function for TaylorModel, not implemented!
template <class T> TaylorModel<T> sinh(const TaylorModel<T>&) {throw("sinh(TaylorModel) not implemented!");}
//! Free cosh() function for TaylorModel, not implemented!
template <class T> TaylorModel<T> cosh(const TaylorModel<T>&) {throw("cosh(TaylorModel) not implemented!");}

namespace std{
//! Multiply complex TaylorModel of type double with rhs TaylorModel of type double.
template <> template <> std::complex<TaylorModel<double> > &std::complex<TaylorModel<double> >::operator*=(const std::complex<TaylorModel<double> > &rhs);
//! Multiply complex TaylorModel of type Iv with rhs TaylorModel of type Iv.
template <> template <> std::complex<TaylorModel<Iv> > &std::complex<TaylorModel<Iv> >::operator*=(const std::complex<TaylorModel<Iv> > &rhs);
}
//! Free function to allow multiplying two const complex TaylorModels
template <class T> std::complex<TaylorModel<T> > operator*(const std::complex<TaylorModel<T> > &t1,const std::complex<TaylorModel<T> > &t2);

#endif //TAYLORMODEL_H
