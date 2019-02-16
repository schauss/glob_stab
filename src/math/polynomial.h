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

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <boost/operators.hpp>

#include "monomial.h"

template <class T>
class Polynomial;

//! Stream operator for debug-output
template <class T>
std::ostream & operator<<(std::ostream &os, const Polynomial<T>& p);

/*!
 * \brief The Polynomial class implements multivariate polynomials polynomials
 *
 * Arithmetic operations as well as functions to calculate an outer
 * approximation of the polynomial using interval arithmetics or Bernstein
 * polynomials are provided.
 *
 * A Polynomial is represented as a set of Monomials and all arithmetic
 * operations are broken down to arithmetic operations on these Monomials.
 */
template <class T=double>
class Polynomial
        : boost::additive< Polynomial<T>
        , boost::multipliable< Polynomial<T>
        , boost::additive2< Polynomial<T>, double
        , boost::multipliable2< Polynomial<T>, double
        > > > >
{
public:
    //! Do Bernstein bound in bound() by default?
    static void doBernsteinBound(bool dobb) {do_bernstein_bound = dobb;}
    //! Do Bernstein bound in bound() by default?
    static bool doBernsteinBound() {return do_bernstein_bound;}

    //! Construct empty Polynomial
    Polynomial();
    //! Construct Polynomial consisting of one variable named var
    Polynomial(const std::string& var);
    //! Construct Polynomial consisting of one variable named var
    Polynomial(const char var[]);
    //! Construct Polynomial consisting of numeric value val
    Polynomial(T val);
    //! Destructor
    ~Polynomial();

    //! Get set of Monomials representing Polynomial
    const std::set<Monomial<T> > &val() const {return monomial;}

    //! Return constant part of Polynomial
    T constant() const;
    /*!
     * \brief Calculate and return outer approximation of Polynomial
     * \param domain Range of all parameters (one interval, i.e., all the same)
     * \param do_bernstein Use Bernstein transformation to approximate bounds.
     * The default value is set using doBernsteinBound(bool).
     * \return Interval representing outer approximation of Polynomial
     */
    const Iv &bound(const Iv &domain, bool do_bernstein = do_bernstein_bound) const;
    /*!
     * \brief Limit total degree of Polynomial and return outer approximation of
     * sum of monomials above total degree.
     * \param total_degree Maximum sum of exponents per Monomial
     * \param domain Interval to consider for each parameter when evaluating the
     * outer approximation of monomials above total degree
     * \return Outer approximation of sum of Monomials whith exponents larger
     * than total_degree as interval.
     *
     * This function does not use a Bernstein transformation to determine the
     * outer approximation but instead simply returns the sum of the ranges of
     * each Monomial above total_degree.
     *
     * \todo Determine Polynomial of Monomials of degree > total_degree and then
     * evaluate outer approximation of this using Bernstein polynomials to get a
     * more exact approximation.
     */
    Iv cutoff(size_t total_degree, const Iv &domain);

    /*!
     * \brief Coordinate transformation of Polynomial
     * \param old_dom Old domain of Polynomial
     * \param new_dom New domain of Polynomial
     * \return Transformed Polynomial
     *
     * Evaluating the transformed Polynomial on the new_dom is the same as
     * evaluationg the original Polynomial on the old_dom.
     */
    Polynomial changeDomain(const Iv &old_dom, const Iv &new_dom) const;
    //! Element-wise maximum order of exponents
    std::vector<uint8_t> maxOrders() const;

    //! Add Polynomial rhs to this Polynomial returning my address
    Polynomial &operator+= (const Polynomial &rhs);
    //! Subtract Polynomial rhs from this Polynomial returning my address
    Polynomial &operator-= (const Polynomial &rhs);
    //! Multiply Polynomial rhs with this Polynomial returning my address
    Polynomial &operator*= (const Polynomial &rhs);

    //! Add double-value to this Polynomial returning my address
    Polynomial &operator+= (double rhs);
    //! Subtract double-value from this Polynomial returning my address
    Polynomial &operator-= (double rhs);
    //! Multiply this Polynomial by double-value returning my address
    Polynomial &operator*= (double rhs);

    //! Negation of Polynomial.
    Polynomial operator-() const;
    //! i-th power of Polynomial.
    Polynomial power(unsigned int i) const;

    //! Check whether Polynomial is empty.
    bool isEmpty() const {return monomial.empty();}

private:
    //! Set of Monomials. The sum of all of these Monomials represents the Polynomial
    std::set<Monomial<T> > monomial;
    //! Interval bound of the Polynomial (cache for bound())
    mutable Iv poly_bound;
    //! Domain for which poly_bound is valid (cache for bound())
    mutable Iv bound_domain;
    //! poly_bound determined using Bernstein transform? (cache for bound())
    mutable bool did_bernstein;

    static bool do_bernstein_bound; //! Default value for bound()

    //! Stream operator for debug-output
    friend std::ostream & operator<< <> (std::ostream &os, const Polynomial<T>& p);
};

#endif //POLYNOMIAL_H
