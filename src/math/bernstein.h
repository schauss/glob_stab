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

#ifndef BERNSTEIN_H
#define BERNSTEIN_H

#include "bernstein_exponent.h"
#include "iv.h"
#include "vec.hpp"

#include <iostream>
#include <inttypes.h>


template <class T>
class Polynomial;

/*!
 * \brief The Bernstein class implements Bernstein polynomials
 *
 * Several important algorithms on Bernstein polynomials are provided.
 */
template <class T>
class Bernstein
{
public:
    //! Construct empty Bernstein polynomial.
    Bernstein();
    /*!
     * \brief Construct Bernstein polynomial from Polynomial.
     * \param pol Polynomial on the domain [0, 1]
     * \param degree Degree of Bernstein polynomial, defaults to maximum degree of Polynomial
     */
    Bernstein(const Polynomial<T> &pol, std::vector<uint8_t> degree = std::vector<uint8_t>());
    /*!
     * \brief Construct Bernstein polynomial from Bernstein coefficients and BernsteinExponent.
     * \param pol Bernstein coefficient vector
     * \param exp BernsteinExponent corresponding to pol
     */
    Bernstein(const Vec<T> pol, BernsteinExponent *exp);
    /*!
     * \brief Copy constructor.
     * \param b Bernstein polynomial
     */
    Bernstein(const Bernstein<T> &b);
    //! Destructor
    ~Bernstein();

    /*!
     * \brief Assign one Bernstein polynomial to another one.
     * \param b Bernstein polynomial
     */
    Bernstein<T>& operator=(const Bernstein<T> &b);

    //! Bernstein polynomial is empty?
    bool isEmpty() const {return empty;}
    //! Upper bound of bernstein polynomial is exact?
    bool upperSharp() const;
    //! Lower bound of bernstein polynomial is exact?
    bool lowerSharp() const;

    //! Return outer approximation of the value set of the Bernstein polynomial.
    Iv bound() const;
    //! Return inner approximation of the value set of the Bernstein polynomial.
    Iv innerBound() const;

    //! Return maximum range between Bernstein coefficients along one dimension while all other indices are fixed.
    Vec<double> range() const;

    //! Return estimation of maximum of partial derivative along each dimension.
    Vec<double> dx() const;
    //! Return index of dimension along which hte estimated maximum of the partial derivative dx() is largest.
    int maxDxDir() const;

    /*!
     * \brief Split the Bernstein polynomial along dimension dim at point split_point resulting in two new Bernstein
     * polynomials (i.e., sets of Bernstein coefficients).
     * \param dim Dimension along which to split
     * \param split_point Point at which to split, 0 < split_point < 1 (0: minimum value, 1: the maximum value)
     */
    Vec<Bernstein<T> > splitAt(size_t dim, double split_point = 0.5) const;

    //! Return maximum number of simultaneously instantiated non-empty Bernstein polynomials.
    static size_t maxCount() {return max_count;}

    //! Return Bernstein coefficients as vector.
    const Vec<T>& coeffs() const {return poly;}
     //! Return Bernstein coefficients of vertices as vector.
    Vec<T> coeffs_vertex() const;

    //! Return order of Bernstein polynomial (multi-dimensional, i.e., each dimension may have a different order).
    Vec<uint8_t> getOrders() const;
    //! Return index of Bernstein coefficient associated with multi-dimensional exponent
    unsigned int getIdx(Vec<uint8_t> &exp) const;

private:
    //! Calculate the Bernstein coefficients from a normalized polynomial.
    void fastCoeff() const;
    //! Calculate outer and inner approximation of the value set of the polynomial.
    void calcBound() const;
    //! Estimate maximum partial derivative along each dimension.
    void calcDerivative() const;
    //! Calculate maximum range of Bernstein coefficients along each dimension
    void calcRange() const;

    /*!
     * \brief Difference between two Bernstein coefficients.
     * \param id1 Index of first Bernstein coefficient
     * \param id2 Index of second Bernstein coefficient
     */
    double deltaMax(size_t id1, size_t id2) const;
    /*!
     * \brief Maximum value of a Bernstein coefficients.
     * \param id Index of the Bernstein coefficient
     */
    double maxVal(size_t id) const;
    /*!
     * \brief Minimum value of a Bernstein coefficients.
     * \param id Index of the Bernstein coefficient
     */
    double minVal(size_t id) const;

    mutable Vec<T> poly; //!< Bernstein coefficients (coeff_done==True) or normalized polynomial coefficients.
    mutable bool coeff_done; //!< True if the values in poly are Bernstein coefficients, False otherwise.

    mutable Iv iv_bound; //!< Outer approximation of the value set of the polynomial
    mutable Iv iv_inner_bound; //!< Inner approximation of the value set of the polynomial
    mutable bool iv_bound_done; //!< Have the inner and outer approximation been computed?

    bool empty; //!< Is the Bernstein polynomial "empty", i.e., does it have no coefficients?

    mutable Vec<double> directional_range; //!< Maximum range of Bernstein coefficients along each dimension.

    mutable Vec<double> derivative; //!< Estimated maximum partial derivative.
    mutable int max_derivative_dir; //!< Direction in which the estimated maximum partial derivative is largest.

    BernsteinExponent *ex; //!< BernsteinExponent of the Bernstein polynomial

    static size_t count; //!< Current number of instantiated non-empty Bernstein polynomials.
    static size_t max_count; //!< Maximum number of simultaneously instantiated non-empty Bernstein polynomials.

    //! Stream operator for debug-output
    template <class U>
    friend std::ostream & operator<<(std::ostream &os, const Bernstein<U>& b);
};

//! Stream operator for debug-output
template <class T>
std::ostream & operator<<(std::ostream &os, const Bernstein<T>& b);

#endif // BERNSTEIN_H
