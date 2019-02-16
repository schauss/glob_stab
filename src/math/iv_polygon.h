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

#ifndef IV_POLYGON_H
#define IV_POLYGON_H

#include <iostream>

#include "polygon.h"
#include "vec.hpp"
#include "iv.h"

/*!
 * \brief The IvPolygon class actually represents a number of 2d-interval boxes.
 *
 * Much of the functionality of the class Polygon is not implemented. This is
 * the reason why this class is not a specialization of the Polygon-class.
 *
 * The main functionality this class offers is:
 *
 * * Construction from vector of \a x and \a y intervals (Iv)
 * * mid(): return Polygon of mid-points of interval boxes.
 * * convexHull(): return the convex hull of all interval-boxes as Polygon.
 *
 * \todo The name is confusing as we would expect this class to represent
 * arbitrary polygons and allow similar operations as Polygon. So we should
 * either change the name or offer the same functionality as Polygon and make
 * this a specialization of the Polygon template class.
 */
class IvPolygon
{
public:
    //! Construct empty IvPolygon.
    IvPolygon();
    //! Construct IvPolygon with interval boxes with given \a x and \a y values.
    IvPolygon(const Vec<Iv> &x, const Vec<Iv> &y);
    //! Destructor.
    ~IvPolygon();

    //! Return number of interval-boxes.
    size_t size() const {return x.size();}
    //! Number of interval-boxes is empty?
    bool empty() {return x.empty();}

    //! Return the convex hull of all interval-boxes as Polygon.
    Polygon convexHull() const;
    //! Return Polygon of mid-points of interval boxes.
    Polygon mid() const;

    //! Append the interval-box \a p (a vector of Iv of length two).
    void push_back(const Vec<Iv> p);
    //! Clear all interval-boxes
    void clear();

    //! Return interval-box \a i (as vector of Iv of length two).
    Vec<Iv> operator[](size_t i) const;

private:
    //! Vector of intervals in x-direction.
    Vec<Iv> x;
    //! Vector of intervals in y-direction
    Vec<Iv> y;

    //! Stream operator for debug-output
    friend std::ostream & operator<<(std::ostream &os, const IvPolygon& p);
};

//! Stream operator for debug-output
std::ostream & operator<<(std::ostream &os, const IvPolygon& p);

#endif //IV_POLYGON_H
