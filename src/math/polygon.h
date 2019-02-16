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

#ifndef POLYGON_H
#define POLYGON_H

#include <iostream>
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#pragma GCC diagnostic pop
#include <CGAL/Polygon_2.h>

#include "vec.hpp"
#include "iv.h"

/*!
 * \brief The Polygon class represents a polygon in the two-dimensional plane.
 *
 * It is used to check whether zero (or an interval) is contained in the convex
 * hull of a set of points. The class makes use of CGAL to perform these checks
 * efficiently and exactly, i.e., an underlying multi-precision numeric type is
 * used within CGAL to determine whether a point is exactly on a polygon or
 * whether it is included or excluded from the polygon.
 *
 * Note that some of the used CGAL-functionality is licensed under the LGPL
 * while other parts are licensed under the GPL.
 *
 * \todo Polygon::Kernel as well as Polygon::Exact_Kernel are defined
 * identically, so the conversions we use in some places don't make any sense
 * (although they don't do any harm) and we could simplify things a little by
 * only using Polygon::Exact_point_2 and Polygon::Exact_Polygon_2.
 */
class Polygon
{
public:
    //typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
    typedef CGAL::Cartesian_converter<Kernel,Exact_Kernel> Inexact_to_Exact;
    typedef CGAL::Cartesian_converter<Exact_Kernel,Kernel> Exact_to_Inexact;

    typedef Kernel::Point_2 Point_2;
    typedef Exact_Kernel::Point_2 Exact_Point_2;
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Polygon_2<Exact_Kernel> Exact_Polygon_2;

    //! Inclusion
    enum include_t {NOT_INCLUDED=0, PARTIALLY_INCLUDED, COMPLETELY_INCLUDED};

    //! Construct an empty polygon.
    Polygon();
    //! Construct a polygon from vectors of x-values and y-values.
    Polygon(const Vec<double> &x, const Vec<double> &y);
    //! Construct a polygon from a vector of 2d-points.
    Polygon(const Vec<Point_2> &p);
    //! Construct a polygon from a Polygon::Exact_Polygon_2 (CGAL-type).
    Polygon(const Exact_Polygon_2 &p);
    //! Destructor.
    ~Polygon();

    //! Cast to Polygon::Exact_Polygon_2 which is a CGAL-type.
    operator Exact_Polygon_2() const;

    //! Add a 2d-point to a Polygon
    const Polygon &push_back(const Point_2 &p);
    //! Return a 2d-point of a polygon
    const Point_2 &operator[](size_t i) const {return poly[i];}
    //! Clear the polygon
    void clear() {poly.clear();}
    //! Number of points in polygon
    size_t size() const {return poly.size();}

    //! Determine the convex hull of the polygon and return it as Polygon
    Polygon convexHull() const;
    //! Return the intersection with another Polygon as Polygon.
    Polygon intersect(const Polygon &p) const;
    //! Polygon is empty, i.e. inner area of polygon is zero.
    bool empty() const;

    //! Check whether x/y is included in the polygon.
    template <class T>
    include_t includes(const T &x, const T &y) const;
    //! Check whether zero is included in the polygon.
    bool includesZero() const {return includes(0.0,0.0);}
    /*!
     * \brief If set to true, the polygon boundary is considered as included.
     *
     * See Polygon::touches_includes for more details.
     */
    void touchesIncludes(bool val) {touches_includes = val;}
    /*!
     * \brief If set to true, the polygon boundary is considered as excluded.
     *
     * See Polygon::touches_excludes for more details.
     */
    void touchesExcludes(bool val) {touches_excludes = val;}

private:
    /*!
     * \brief Considers polygon to be simple and checks whether the value x/y is
     * included in the polygon.
     *
     * See simplify() for more information on simple Polygons.
     */
    template <class T>
    include_t simpleIncludes(const T &x, const T &y) const;
    /*!
     * \brief Considers polygon to be simple and checks whether the simple
     * Polygon \a p is included in the polygon.
     * \param p simple Polygon
     *
     * See simplify() for more information on simple polygons.
     */
    include_t simpleIncludes(const Polygon &p) const;

    //! Check whether the bounding box of the polygon contains x/y.
    template <class T>
    bool boxIncludes(const T &x, const T &y) const;

    /*!
     * \brief Simplify a non-simple polygon.
     *
     * A simple polygon is a polygon where no line-segments intersect.
     * Non-simple polygons can generally be represented by a number of simple
     * polygons.
     */
    void simplify() const;
    Vec<Polygon> makeTriangles(int idx) const;

    //! Stores the polygon as CGAL-type
    Polygon_2 poly;
    //! Cache of simple polygons. Created by simplify() if polygon is not simple.
    mutable Vec<Polygon> simple;
    //! Is this a simple polygon? Set to true by simplify() if polygon is simple.
    mutable bool is_simple;

    /*!
     * \brief If true, a point on the boundary is considered as being included
     * in the Polygon.
     *
     * touches_excludes can be set independently!
     */
    bool touches_includes;
    /*!
     * \brief If true, a point on the boundary is considered as being excluded
     * in the Polygon.
     *
     * touches_includes can be set independently!
     */
    bool touches_excludes;

    //! Stream operator for debug-output
    friend std::ostream & operator<<(std::ostream &os, const Polygon& p);
};

//! Stream operator for debug-output
std::ostream & operator<<(std::ostream &os, const Polygon& p);

#endif //POLYGON_H
