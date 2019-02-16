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

#ifndef BERNSTEIN2D_H
#define BERNSTEIN2D_H

#include "vec.hpp"
#include "bernstein.h"
#include "bernstein_2d_patch.h"
#include "polygon.h"
#include "iv_polygon.h"
#include "iv.h"
#include "taylor_model.h"

/*!
 * \brief The Bernstein2d class implements two-dimensional multivariate Bernstein polynomials.
 *
 * A two-dimensional multivariate Bernstein polynomial is a function which maps a number of variables to two values
 * using Bernstein polynomials. This can, e.g., be used to represent complex-values polynomials which is necessary for a
 * stability analysis using the value-set approach / zero-exclusion principle.
 *
 * Different operations can be performed on the 2d Bernstein polynomial, e.g., split the domain, calculate the estimated
 * maximum partial derivative and its direction, and evaluate zero inclusion/exclusion as well as inclusion/exclusion
 * of an interval rest in the inner and outer approximation of the value set of the 2D Bernstein polynomial.
 *
 * \todo When calling functions like restExcluded, zeroIncluded, etc. we should check whether the evaluate was called
 * and call it otherwise!
 */
template <class T=double>
class Bernstein2d
{
public:
    //! Construct empty 2D Bernstein polynomial.
    Bernstein2d();
    /*!
     * \brief Construct 2D Bernstein polynomial from a TaylorModels.
     * \param tm_res Vector of TaylorModels (length two required)
     */
    Bernstein2d(const Vec<TaylorModel<T> > &tm_res);
    /*!
     * \brief Construct 2D Bernstein polynomial from Bernstein polynomials and interval remainders.
     * \param b Vector of Bernstein polynomials (length two required)
     * \param rest_iv Vector of interval remainders (length two required)
     */
    Bernstein2d(const Vec<Bernstein<T> > &b, const Vec<Iv> rest_iv);
    //! Destructor
    ~Bernstein2d();

    //! Evaluate inner and outer approximations of 2D Bernstein polynomial and check for zero exclusion and inclusion
    void evaluate();

    //! Store the convex hull to the Matlab m-file b_coeffs.m
    void storeConvexHull() const;

    /*!
     * \brief Split the 2d Bernstein polynomial along dimension dim at point split_point resulting in two new 2d
     * Bernstein polynomials. This uses the Bernstein::splitAt() on both Bernstein polynomials and returns two 2d
     * Bernstein polynomials.
     * \param dim Dimension along which to split
     * \param split_point Point at which to split, 0 < split_point < 1 (0: minimum value, 1: the maximum value)
     */
    Vec<Bernstein2d> splitAt(size_t dim, double split_point = 0.5) const;

    //! Calculate and return the sum of square estimated maximum partial derivatives in each direction.
    Vec<double> dx() const;
    //! Return the dimension of dx() with maximum value
    size_t maxDxDir() const;

    bool restIncluded() const;            //!< Negative interval remainder included, i.e., zero included in function
    bool restExcluded() const;            //!< Negative interval remainder excluded, i.e., zero excluded from function
    bool restPartiallyIncluded() const;   //!< Negative interval remainder partly included
    bool restPartiallyExcluded() const;   //!< Negative interval remainder partly excluded
    bool zeroIncluded() const;            //!< Zero included, i.e., zero included in polynomial part
    bool zeroExcluded() const;            //!< Zero excluded, i.e., zero exluded from polynomial part
    bool unknown() const;                 //!< Inclusion and exclusion is not known (i.e., no sufficient condition holds)
    bool isEmpty() const;                 //!< 2d Bernstein polynomial is empty

    bool boundaryIsZero() const;          //!< The boundary of the convex hull passes through zero exactly

private:
    void setupPoints() const;             //!< Setup polygon from Bernstein coefficients
    void computeConvexHull() const;       //!< Evaluate convex hull of all Bernstein coefficients
    void checkPatches() const;            //!< Check zero inclusion
    void checkPatchesIv() const;          //!< Check zero inclusion for Interval Bernstein coefficients

    Vec<Bernstein<T> > poly;              //!< Vector of Bernstein polynomials (length two)
    Vec<Iv> rest;                         //!< Vector of Interval remainders (length two)
    mutable Polygon points;               //!< Polygon containing the 2d points of all Bernstein coefficients
    mutable IvPolygon iv_points;          //!< Polygon containing the 2d intervals of all Bernstein coefficients
    mutable Polygon hull;                 //!< Polygon representing the covex hull of all Bernstein coefficients

    mutable bool rest_included;           //!< Negative interval remainder included, i.e., zero included in function
    mutable bool rest_excluded;           //!< Negative interval remainder excluded, i.e., zero excluded from function
    mutable bool rest_partially_included; //!< Negative interval remainder partly included
    mutable bool rest_partially_excluded; //!< Negative interval remainder partly excluded
    mutable bool zero_included;           //!< Zero included, i.e., zero included in polynomial part
    mutable bool zero_excluded;           //!< Zero excluded, i.e., zero exluded from polynomial part
    bool empty;                           //!< 2d Bernstein polynomial is empt

    //! Stream operator for debug-output
    template <class U>
    friend std::ostream & operator<<(std::ostream &os, const Bernstein2d<U>& b2d);
};

//! Stream operator for debug-output
template <class T>
std::ostream & operator<<(std::ostream &os, const Bernstein2d<T>& b2d);

#endif // BERNSTEIN2D_H
