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

#include "polygon.h"

#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/intersections.h>
#include <CGAL/Boolean_set_operations_2.h>
#pragma GCC diagnostic pop

Polygon::Polygon()
    : is_simple(false), touches_includes(false), touches_excludes(false)
{

}

Polygon::Polygon(const Vec<double> &x, const Vec<double> &y)
    : is_simple(false), touches_includes(false), touches_excludes(false)
{
    if (x.size() != y.size())
        throw("Polygon::Polygon(x,y) - different sizes for x and y not allowed!");

    //poly.container().reserve(x.size());

    for (size_t i=0; i<x.size(); i++)
        poly.push_back(Point_2(x[i],y[i]));
}

Polygon::Polygon(const Vec<Point_2> &p)
    : is_simple(false), touches_includes(false), touches_excludes(false)
{
    poly = Polygon_2(p.begin(), p.end());
}

Polygon::Polygon(const Exact_Polygon_2 &p)
    : is_simple(false), touches_includes(false), touches_excludes(false)
{
    Exact_to_Inexact etoi;

    for (size_t i=0; i<p.size(); i++)
        push_back(etoi(p[i]));
}

Polygon::~Polygon()
{

}

Polygon::operator Exact_Polygon_2() const
{
    Exact_Polygon_2 p;

    Inexact_to_Exact itoe;

    for (size_t i=0; i<poly.size(); i++)
        p.push_back(itoe(poly[i]));

    return p;
}

const Polygon& Polygon::push_back(const Point_2 &p)
{
    poly.push_back(p);
    is_simple = false;
    simple.clear();
    return *this;
}

/*!
 * To determine the convex hull a Graham-Andrew's Scan implemented in CGAL is
 * used. See CGAL::ch_graham_andrew for more details.
 * \todo This is a simple polygon! Should set this, then it can be checked in
 * simpleIncludes()!
 */
Polygon Polygon::convexHull() const
{
    Vec<Point_2> hull = Vec<Point_2>(poly.size(),Point_2(0.0,0.0));
    Vec<Point_2>::iterator hull_end = CGAL::ch_graham_andrew(poly.container().begin(),poly.container().end(),hull.begin());
    hull.resize(hull_end-hull.begin());
    return hull;
}

/*!
 * The resulting polygon may contain holes, i.e., it is not necessarily simple!
 */
Polygon Polygon::intersect(const Polygon &p) const
{
    if ( (poly.size() < 3) || (p.poly.size() < 3) )
        return Polygon();

    Exact_Polygon_2 a = *this;
    if (a.is_clockwise_oriented())
        a.reverse_orientation();

    Exact_Polygon_2 b = p;
    if (b.is_clockwise_oriented())
        b.reverse_orientation();

    Vec<CGAL::Polygon_with_holes_2<Exact_Kernel> > isec;
    try{
    intersection(a,b,std::back_inserter(isec));
    }
    catch (std::exception &err) {
        std::cerr << "std::exception in CGAL-intersection: " << std::endl << err.what() << std::endl;
        throw;
    }

    if (isec.empty()) {
        return Polygon();
    }
    else {
        if (isec.size() > 1)
            throw("Polygon::intersect results in two disjoint polygonials!");
        if (isec[0].has_holes())
            throw("Polygon::intersect - has holes!");
        if (isec[0].is_unbounded())
            throw("Polygon::intersect - is unbounded!");

        Exact_Polygon_2 ret_exact = isec[0].outer_boundary();

        return Polygon(ret_exact);
    }
}

bool Polygon::empty() const
{
    if (poly.is_empty())
        return true;

    simplify();

    if (is_simple) {
        return (poly.area() == 0.0);
    } else {
        return ( apply(simple,&Polygon::empty).all() );
    }
}

/*!
 * \param x,y Intervals (Iv) or double
 * \return Inclusion as Polygon::include_t
 *
 * Simplify polygon and then either check polygon (if it is simple) or all
 * simple sub-polygons for inclusion of x/y.
 *
 * The return value is as evaluated with this prescedence:
 *
 * * x/y is COMPLETELY_INCLUDED in one simple polygon => COMPLETELY_INCLUDED
 * * x/y is PARTIALLY_INCLUDED  in one simple polygon => PARTIALLY_INCLUDED
 * * x/y is NOT_INCLUDED in any simple polygon => NOT_INCLUDED
 *
 * See the functions which are called for the simple polygons for a detailed
 * description of the return values:
 *
 * * \a x, \a y doubles: simpleIncludes(const T &x, const T &y) const
 * * \a x, \a y intervals (Iv): simpleIncludes(const Iv &x, const Iv &y) const
 */
template <class T>
Polygon::include_t Polygon::includes(const T &x, const T &y) const
{
    // first perform a simple check (box-inclusion)
    if (!boxIncludes(x,y))
        return NOT_INCLUDED;

    simplify();

    if (is_simple) {
        return simpleIncludes(x,y);
    } else {
        // check all polygons in simple
        //  - if any completely includes (x,y) return completely
        //  - else if any partially includes (x,y) return partially
        //  - else return not
        include_t ret = NOT_INCLUDED;
        for (size_t i = 0; i<simple.size(); i++) {
            include_t inc = simple[i].simpleIncludes(x,y);
            if (inc == COMPLETELY_INCLUDED)
                return COMPLETELY_INCLUDED;
            else if (inc == PARTIALLY_INCLUDED)
                ret = PARTIALLY_INCLUDED;
        }
        return ret;
    }
}

/*!
 * \relates Polygon
 * \brief Explicit instantiation of includes(const T &x, const T &y) const for
 * type double.
 */
template Polygon::include_t Polygon::includes<double>(const double &x, const double &y) const;
/*!
 * \relates Polygon
 * \brief Explicit instantiation of includes(const T &x, const T &y) const for
 * interval-type Iv.
 */
template Polygon::include_t Polygon::includes<Iv>(const Iv &x, const Iv &y) const;

/*!
 * The return values are:
 *
 * * NOT_INCLUDED if \a x/\a y is completely excluded from the polygon.
 * * COMPLETELY_INCLUDED if \a x/\a y is included in the polygon.
 * * PARTIALLY_INCLUDED if \a x/\a y is is on the boundary of the polygon.
 *
 * See simpleIncludes(const Iv&, const Iv&) const for the specialization for
 * interval-type Iv.
 */
template <class T>
Polygon::include_t Polygon::simpleIncludes(const T &x, const T &y) const
{
    // removed check as this is quite expensive, and simpleIncludes is only called from includes (after simplify)
    //if (!poly.is_simple())
    //    throw("Polygon::simpleIncludes(double,double) called for non-simple polygon!");

    switch( CGAL::bounded_side_2(poly.vertices_begin(), poly.vertices_end(), Point_2(x,y)) ) {
    case CGAL::ON_BOUNDED_SIDE :
        //std::cout << " Polygon includes (" << x << ", " << y << ")" << std::endl;
        return COMPLETELY_INCLUDED;
    case CGAL::ON_BOUNDARY:
        //std::cout << " Polygon touches (" << x << ", " << y << ")" << std::endl;
        return PARTIALLY_INCLUDED;
    case CGAL::ON_UNBOUNDED_SIDE:
    default:
        return NOT_INCLUDED;
    }
}

/*!
 * Specialization of simpleIncludes(const T&, const T&) const for type Iv.
 *
 * Evalation is a little more complex in this case. The return values are:
 *
 * * NOT_INCLUDED if the interval-box \a x/\a y is completely excluded from the
 *   polygon.
 * * COMPLETELY_INCLUDED if the interval-box \a x/\a y is completely included in
 *   the polygon.
 * * PARTIALLY_INCLUDED otherwise, i.e., if the interal-box intersects the
 *   polygon boundary.
 *
 * The special cases where the interval-box is included but touches the polygon
 * boundary and where the interval-box is excluded but touches the polygon
 * boundary are handled depending on the values of Polygon::touches_includes and
 * Polygon::touches_excludes:
 *
 * * interval-box included in polygon but touching boundary
 *   * Polygon::touches_includes==true => COMPLETELY_INCLUDED
 *   * Polygon::touches_includes==false => PARTIALLY_INCLUDED
 * * interval-box excluded from polygon but touching boundary
 *   * Polygon::touches_excludes==true => NOT_INCLUDED
 *   * Polygon::touches_excludes==false => PARTIALLY_INCLUDED
 *
 * This function actually handles the degenerate case where the interval-box
 * consists of a line and calls simpleIncludes(const T&, const T&) const in case
 * the interval-box is actually a point and simpleIncludes(const Polygon&) const
 * otherwise.
 */
template <>
Polygon::include_t Polygon::simpleIncludes(const Iv &x, const Iv &y) const
{
    // removed check as this is quite expensive, and simpleIncludes is only called from includes (after simplify)
    //if (!poly.is_simple())
    //    throw("Polygon::simpleIncludes(Iv,Iv) called for non-simple polygon!");

    if ( (x.width() == 0.0) && (y.width() == 0.0) ){
        return simpleIncludes(x.lower(), y.lower());
    } else if ( (x.width() == 0.0) || (y.width() == 0.0) ) {
        // degenerated polygon, i.e. a line - first check both endpoints
        include_t a = simpleIncludes(x.lower(),y.lower());
        include_t b = simpleIncludes(x.upper(),y.upper());
        bool one_included = ( (a==COMPLETELY_INCLUDED) || (b==COMPLETELY_INCLUDED) );
        bool one_excluded = ( (a==NOT_INCLUDED) || (b==NOT_INCLUDED) );
        bool one_touches  = ( (a==PARTIALLY_INCLUDED) || (b==PARTIALLY_INCLUDED) );

        //if only one endpoint is included, the line must intersect the polynomial-boundary!
        if ( (!one_included && !one_excluded) || ( one_included && one_excluded ) || (!touches_includes && one_touches && one_included) || (!touches_excludes && one_touches && one_excluded) )
            return PARTIALLY_INCLUDED;
        else {
            // otherwise, we must create a segment connecting both edges and check for intersection with one of the edges of the polynomial
            Inexact_to_Exact itoe;
            Exact_Kernel::Segment_2 iv_line(itoe(Point_2(x.lower(),y.lower())),itoe(Point_2(x.upper(),y.upper())));
            for (size_t i=0; i<poly.size(); i++) {
                CGAL::Object obj = intersection(itoe(poly.edge(i)),iv_line);
                if ( const Exact_Point_2 *p_int = CGAL::object_cast<Exact_Point_2>(&obj) ) {
                    if ( (!touches_includes && !touches_excludes) || ( (*p_int != iv_line.point(0)) && (*p_int != iv_line.point(1)) ) )
                        return PARTIALLY_INCLUDED;
                }
            }

            // either both endpoints are included or none, and the polygon is not intersected
            if (one_included)
                return COMPLETELY_INCLUDED;
            else if (!one_included)
                return NOT_INCLUDED;
            else
                throw("Polygon::simpleIncludes(Iv,Iv) - Cannot determine inclusion!");
        }
    } else {
        Polygon p;
        p.push_back(Point_2(x.lower(),y.lower()));
        p.push_back(Point_2(x.upper(),y.lower()));
        p.push_back(Point_2(x.upper(),y.upper()));
        p.push_back(Point_2(x.lower(),y.upper()));

        return simpleIncludes(p);
    }
}

/*!
 * \relates Polygon
 * \brief Explicit instantiation of simpleIncludes(const T &x, const T &y) const
 * for type double.
 */
template Polygon::include_t Polygon::simpleIncludes<double>(const double &x, const double &y) const;
/*!
 * \relates Polygon
 * \brief Explicit instantiation of simpleIncludes(const T &x, const T &y) const
 * for type Polygon::Kernel::FT.
 */
template Polygon::include_t Polygon::simpleIncludes<Polygon::Kernel::FT>(const Polygon::Kernel::FT &x, const Polygon::Kernel::FT &y) const;
/*!
 * \relates Polygon
 * \brief Explicit instantiation of simpleIncludes(const T &x, const T &y) const
 * for type Iv.
 */
template Polygon::include_t Polygon::simpleIncludes<Iv>(const Iv &x, const Iv &y) const;

/*!
 * The return values are:
 *
 * * NOT_INCLUDED if the \a p is completely excluded from the polygon.
 * * COMPLETELY_INCLUDED if \a p is completely included in the polygon.
 * * PARTIALLY_INCLUDED otherwise, i.e., if \a p intersects the polygon
 * boundary.
 *
 * The special cases where \a p is included but touches the polygon boundary and
 * where \a p is excluded but touches the polygon boundary are handled depending
 * on the values of Polygon::touches_includes and Polygon::touches_excludes:
 *
 * * \a p included in polygon but touching boundary
 *   * Polygon::touches_includes==true => COMPLETELY_INCLUDED
 *   * Polygon::touches_includes==false => PARTIALLY_INCLUDED
 * * \a p excluded from polygon but touching boundary
 *   * Polygon::touches_excludes==true => NOT_INCLUDED
 *   * Polygon::touches_excludes==false => PARTIALLY_INCLUDED
 */
Polygon::include_t Polygon::simpleIncludes(const Polygon &p) const
{
    // removed check as this is quite expensive, and simpleIncludes is only called from includes (after simplify)
    //if (!poly.is_simple())
    //    throw("Polygon::simpleIncludes(Polygon) called for non-simple polygon!");

    //if (!p.poly.is_simple())
    //    throw("Polygon::simpleIncludes(Polygon p) called for non-simple polygon p!");

    /* for each vertex in p, check if it is included in this:
       if one point is included (->any_included) and one point is not (->!all_included)
       the polygon must be partially included */
    bool one_included = false;
    bool one_excluded = false;
    bool one_touches = false;
    for (size_t i=0; i<p.poly.size(); i++) {
        include_t inc = simpleIncludes(p.poly[i].x(),p.poly[i].y());
        if (inc == NOT_INCLUDED) {
            one_excluded = true;
        } else if (inc == COMPLETELY_INCLUDED) {
            one_included = true;
        } else {
            one_touches = true;
        }
    }

    if (!one_included && !one_excluded)
        std::cout << "ALL TOUCHING!" << std::endl;

    if ( (!one_included && !one_excluded) || ( one_included && one_excluded ) || (!touches_includes && one_touches && one_included) || (!touches_excludes && one_touches && one_excluded) )
        return PARTIALLY_INCLUDED;

    /* if we get here, then either all points are included or all points are not included
       therefore, we check intersection of each segment of p with each segment of this
       - if there is an intersection, p is partially included
       - otherwise, p is either completely included or not included at all */

    // first convert polygon p to exact polygon as we need exact segments (for both polygons) and otherwise
    // we would perform conversion up to poly.size() times for p.poly()
    Exact_Polygon_2 b = p;

    Inexact_to_Exact itoe;
    for (size_t i=0; i<poly.size(); i++) {
        // convert segments of poly on the fly, as this could save conversions
        Exact_Kernel::Segment_2 line1(itoe(poly.edge(i)));
        for (size_t j=0; j<b.size(); j++) {
            CGAL::Object obj = intersection(b.edge(j),line1);
            if ( const Exact_Point_2 *p_int = CGAL::object_cast<Exact_Point_2>(&obj) ) {
                if ( (!touches_includes && !touches_excludes) || ( (*p_int != b[j]) && (*p_int != b[j==b.size()-1?0:j+1]) ) ) {
                    return PARTIALLY_INCLUDED;
                }
            }
        }
    }

    /* if we get here, then there is no intersection of segments, i.e.,
       either the polygon is completely included or not at all */
    if ( (!one_touches || touches_excludes) && !one_included) {
        if (p.simpleIncludes(poly[0].x(),poly[0].y())) // check if this is completely in p!
            return PARTIALLY_INCLUDED;
        else
            return NOT_INCLUDED;
    }
    else if ( (!one_touches || touches_includes) && !one_excluded)
        return COMPLETELY_INCLUDED;
    else
        throw("Polygon::simpleIncludes(Polygon) - Cannot determine inclusion!");
}

/*!
 * This is a necessary condition for the polygon to include x/y.
 *
 * \bug This may not be correct as it only returns true if the bounding box is
 * larger than x/y! What about equality?
 */
template <class T>
bool Polygon::boxIncludes(const T &x, const T &y) const
{
    // in x-direction
    bool max_larger = true;
    bool min_smaller = true;
    for (size_t i=0; i<poly.size(); i++) {
        if (!max_larger && (poly[i].x() > x) ) {
            max_larger = true;
            if (min_smaller)
                break;
        }
        if (!min_smaller && (poly[i].x() < x) ) {
            min_smaller = true;
            if (max_larger)
                break;
        }
    }
    if (!max_larger || !min_smaller)
        return false;

    // in y-direction
    max_larger = true;
    min_smaller = true;
    for (size_t i=0; i<poly.size(); i++) {
        if (!max_larger && (poly[i].y() > y) ) {
            max_larger = true;
            if (min_smaller)
                break;
        }
        if (!min_smaller && (poly[i].y() < y) ) {
            min_smaller = true;
            if (max_larger)
                break;
        }
    }
    if (!max_larger || !min_smaller)
        return false;

    return true;
}

/*!
 * Specialization of boxIncludes(const T&,const T&) const for type Iv.
 *
 * \todo Is this really necessary? Comparison with an Iv should return the
 * correct result, shouldn't it?
 */
template <>
bool Polygon::boxIncludes(const Iv &x, const Iv &y) const
{
    // in x-direction
    bool max_larger = true;
    bool min_smaller = true;
    for (size_t i=0; i<poly.size(); i++) {
        if (!max_larger && (poly[i].x() > x.upper()) ) {
            max_larger = true;
            if (min_smaller)
                break;
        }
        if (!min_smaller && (poly[i].x() < x.lower()) ) {
            min_smaller = true;
            if (max_larger)
                break;
        }
    }
    if (!max_larger || !min_smaller)
        return false;

    // in y-direction
    max_larger = true;
    min_smaller = true;
    for (size_t i=0; i<poly.size(); i++) {
        if (!max_larger && (poly[i].y() > y.upper()) ) {
            max_larger = true;
            if (min_smaller)
                break;
        }
        if (!min_smaller && (poly[i].y() < y.lower()) ) {
            min_smaller = true;
            if (max_larger)
                break;
        }
    }
    if (!max_larger || !min_smaller)
        return false;

    return true;
}

/*!
 * This function first checks whether the polygon is simple and indicates this
 * in the variable Polygon::is_simple.
 *
 * If this is not the case the polygon is simplified and the resulting simple
 * polygons are collected in Polygon::simple. Note that this is not a general
 * implementation but instead constrained to two special cases
 * * convex Polygons: if a polygon is convex but non-simple this actually means
 *   it is a simple polygon where some vertex consists of multiple points.
 * * Co-linear polygon, i.e., all vertices are on one line.
 * * Polygon with 4 distinct vertices where two line-segments intersect. This is
 *   split into two triangles using makeTriangles().
 *
 * This actually covers all possible cases for a polygon with 4 vertices which
 * is the use case of this function.
 */
void Polygon::simplify() const
{
    if (is_simple || !simple.empty())
        return;

    try{
    if ( poly.is_simple()) {
        // simple polygon -> return this
        is_simple = true;

        // removed orientation-fixing, as only necessary when intersecting polygons and done there anyway (and because it is slow!)
        //if (simple.back().poly.is_clockwise_oriented())
        //    simple.back().poly.reverse_orientation();
    } else
        if (poly.is_convex()) {
        // convex polygon -> return convex hull which is simple!
        // this means we either have one triangle or all points are colinear
        simple.push_back(convexHull());
    } else {
        // this means we have two triangles -> split

        simple = makeTriangles(0);
        if (simple.empty()) {
            simple = makeTriangles(1);
            if (simple.empty()) {
                throw("Polygon::simplify - no triangles found when intersecting!");
            }
        }
    }
    } catch (std::exception &err) {
        std::cerr << "std::exception in Polygon::simplify()" << std::endl << err.what() << std::endl;
        throw;
    }
}

/*!
 * \param idx First edge used for intersection
 * \return Vector of triangels, may be empty if requested triangels do not exist
 * or a line if all points are on one line!
 *
 * Intersects the polygon edges idx and idx+2.If these intersect at a point
 * this means the polygon consists of two triangles where each of these
 * triangles consists of one of the other edges connected to the
 * intersection-point. These two triangles are then returned. If the two edges
 * intersect as a segment this means the four points are on one line and the
 * convex hull of the four points is returned as only element of the vector.
 * If the two lines do not intersect an empty vector is returned.
 *
 * \todo Refactor. The check of different cases from Polygon::simplify() could
 * be moved into this function (makeTriangles would not need any arguments).
 */
Vec<Polygon> Polygon::makeTriangles(int idx) const
{
    if(poly.container().size() != 4)
        throw("Polygon::makeTriangle - Polygon does not have 4 vertices!");

    if( (idx < 0) ||(idx > 1) )
        throw("Polygon::makeTriangle - idx out of bounds!");

    Vec<Polygon> triangles;

    try {
    Inexact_to_Exact itoe;
    Exact_Kernel::Segment_2 edge_a = itoe(poly.edge(idx));
    Exact_Kernel::Segment_2 edge_b = itoe(poly.edge(idx+2));
    CGAL::Object obj = intersection(edge_a,edge_b);

    if (const Exact_Point_2 *point = CGAL::object_cast<Exact_Point_2>(&obj)) {
        //std:: cout << "Point intersection between edge " << idx << " and " << (idx+2) << std::endl;

        Exact_to_Inexact etoi;

        Polygon p;
        p.push_back(poly[idx]);
        p.push_back(etoi(*point));
        p.push_back(poly[(idx+3)%4]);

        p.simplify();
        if (!p.is_simple) {
            if (p.simple.size() != 1)
                throw("Polygon::makeTriangle - simplification of triangle 1 failed!");
            p = p.simple[0];
        }
        triangles.push_back(p);

        p.clear();
        p.push_back(poly[idx+1]);
        p.push_back(etoi(*point));
        p.push_back(poly[idx+2]);

        p.simplify();
        if (!p.is_simple) {
            if (p.simple.size() != 1)
                throw("Polygon::makeTriangle - simplification of triangle 2 failed!");
            p = p.simple[0];
        }
        triangles.push_back(p);

    } else if (CGAL::object_cast<Exact_Kernel::Segment_2>(&obj)) {
        // if there is a segment intersection between the two edges, then this must actually be a line, i.e. the convex hull returns the correct result
        triangles.push_back(convexHull());
        //throw("Polygon::makeTriangle - Segment intersection between two edges");
    } else {
        //std:: cout << "No intersection between edge " << idx << " and " << (idx+2) << std::endl;
    }

    } catch (std::exception &err) {
        std::cerr << "std::exception in Polygon::makeTriangles(" << idx << "):" << std::endl << err.what() << std::endl;
        throw;
    }
    return triangles;
}

std::ostream & operator<<(std::ostream &os, const Polygon& p)
{
    os<< Vec<Polygon::Point_2>(p.poly.container()) << std::endl;
    return os;
}
