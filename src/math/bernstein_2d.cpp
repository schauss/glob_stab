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

#include "bernstein_2d.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "helpers.h"

template <class T>
Bernstein2d<T>::Bernstein2d() :
    rest_included(false), rest_excluded(false), rest_partially_included(false), rest_partially_excluded(false), zero_included(false), zero_excluded(false), empty(true)
{

}

template <class T>
Bernstein2d<T>::Bernstein2d(const Vec<TaylorModel<T> > &tm_res) :
    rest_included(false), rest_excluded(false), rest_partially_included(false), rest_partially_excluded(false), zero_included(false), zero_excluded(false), empty(false)
{
    if (tm_res.size() != 2)
        throw("Bernstein2d::Bernstein2d - wrong size of argument!");

    // we must create two bernstein-polynomials with identical degree so that coefficients correspond!
    std::vector<uint8_t> degree = tm_res[0].polynomial_01().maxOrders();
    for (size_t i=0; i<degree.size(); i++) {
        if (degree[i] < tm_res[1].polynomial_01().maxOrders()[i]) {
            degree[i] = tm_res[1].polynomial_01().maxOrders()[i];
        }
    }

    for (size_t i=0; i<tm_res.size(); i++) {
        poly.push_back(Bernstein<T>(tm_res[i].polynomial_01(),degree));
        poly[i].bound();
        rest.push_back(tm_res[i].restBound());
    }
}

/*!
 * BEWARE: This constructor assumes that the two Bernstein polynomials have identical degree but this is no verified!
 */
template <class T>
Bernstein2d<T>::Bernstein2d(const Vec<Bernstein<T> > &b, const Vec<Iv> rest_iv) :
    rest_included(false), rest_excluded(false), rest_partially_included(false), rest_partially_excluded(false), zero_included(false), zero_excluded(false), empty(false)
{
    if ( (b.size() != 2) || (rest_iv.size() != 2) )
        throw("Bernstein2d::Bernstein2d - wrong size of arguments!");

    poly = b;
    rest = rest_iv;
}

template <class T>
Bernstein2d<T>::~Bernstein2d()
{

}

/*!
 * In a first step the convex hull of the Bernstein coefficients is calculated. Then we evaluate whether zero is
 * excluded from the convex hull. If this is true we evaluate whether the negative interval remainder is excluded from
 * the convex hull. For Taylor Models consisting of a Polynomial and interval remainder this implies that zero is
 * excluded from the exact function for which this Taylor Model was evaluated.
 *
 * Next (if zero is not excluded from the function) we evaluate zero inclusion. Therefore, we first check whether the
 * negative interval remainder is included in the convex hull as this is a necessary condition for the following zero
 * inclusion check to succeed. Then, we evaluate zero inclusion (or more specifically inclusion of the negative interval
 * remainder) using the method described in Schauss, Peer, Buss, TAC 2017.
 */
template <class T>
void Bernstein2d<T>::evaluate()
{
    if (isEmpty())
        throw("Bernstein2d::evaluate called on empty Bernstein2d!");

    if (!rest[0].isFinite() || !rest[1].isFinite()) {
        std::cout << "Bernstein2d::evaluate - rest is not finite:" << std::endl << rest << std::endl;
        return;
    }

    computeConvexHull();

    if (!hull.includesZero()) {
        zero_excluded = true;
    }

    // check if inverted interval (-rest) is included, as this means zero is included in polygon + rest
    hull.touchesIncludes(true);
    Polygon::include_t inc;

    if (!points.empty() || !iv_points.empty())
        inc = hull.includes(-rest[0],-rest[1]);
    else
        throw("Bernstein2d::evaluate() - no points setup!");

    if (inc == Polygon::COMPLETELY_INCLUDED) {
        if (!points.empty())
            checkPatches();
        else if (!iv_points.empty())
            checkPatchesIv();
    } else if (inc == Polygon::PARTIALLY_INCLUDED) {
        rest_partially_excluded = true;
    } else if (inc == Polygon::NOT_INCLUDED) {
        rest_excluded = true;
    } else
        throw("Bernstein2d::evaluate - unknown return value!");
}

template <class T>
bool Bernstein2d<T>::restIncluded() const
{
    if (isEmpty())
        throw("Bernstein2d::restIncluded called on empty Bernstein2d!");

    return rest_included;
}

template <class T>
bool Bernstein2d<T>::restExcluded() const
{
    if (isEmpty())
        throw("Bernstein2d::restNotIncluded called on empty Bernstein2d!");

    return rest_excluded;
}

template <class T>
bool Bernstein2d<T>::restPartiallyIncluded() const
{
    if (isEmpty())
        throw("Bernstein2d::restPartiallyIncluded called on empty Bernstein2d!");

    return rest_partially_included;
}

template <class T>
bool Bernstein2d<T>::restPartiallyExcluded() const
{
    if (isEmpty())
        throw("Bernstein2d::restPartiallyExcluded called on empty Bernstein2d!");

    return rest_partially_excluded;
}

template <class T>
bool Bernstein2d<T>::zeroIncluded() const
{
    if (isEmpty())
        throw("Bernstein2d::zeroIncluded called on empty Bernstein2d!");

    return zero_included;
}

template <class T>
bool Bernstein2d<T>::zeroExcluded() const
{
    if (isEmpty())
        throw("Bernstein2d::zeroExcluded called on empty Bernstein2d!");

    return zero_excluded;
}

template <class T>
bool Bernstein2d<T>::unknown() const
{
    if (isEmpty())
        throw("Bernstein2d::isUnknown called on empty Bernstein2d!");

    return ( (!rest_included) && (!rest_excluded) && (!rest_partially_included) && (!rest_partially_excluded) && (!zero_included) );
}

template <class T>
bool Bernstein2d<T>::isEmpty() const
{
    return empty;
}

/*!
 * \todo Currently always returns false, must be reimplemented correctly!
 */
template <class T>
bool Bernstein2d<T>::boundaryIsZero() const
{
    throw("Bernstein2d::boundaryIsZero currently not implemented!");

    if (isEmpty())
        throw("Bernstein2d::boundaryIsZero called on empty Bernstein2d!");

    if (hull.size() == 0) {
        std::cout << "Bernstein2d::boundaryIsZero called for empty hull!" << std::endl;
        return false;
    }

    Vec<double> min_val(2,INFINITY);
    Vec<double> max_val(2,-INFINITY);
    /*
    for (size_t i=0; i<hull.size(); i++) {
        Polygon::Point_2 p = hull[i];
        if (p.x() < min_val[0])
            min_val[0] = p.x();
        if (p.x() > max_val[0])
            max_val[0] = p.x();
        if (p.y() < min_val[1])
            min_val[1] = p.y();
        if (p.y() > max_val[1])
            max_val[1] = p.y();
    }
    */

    if ( (min_val[0] == 0.0) || (min_val[1] == 0.0) || (max_val[0] == 0.0) || (max_val[1] == 0.0) )
        return true;
    else
        return false;
}

template <class T>
Vec<Bernstein2d<T> > Bernstein2d<T>::splitAt(size_t dim, double split_point) const
{
    if (isEmpty())
        throw("Bernstein2d::splitAt() called on empty Bernstein2d!");

    // 2x2 Bernstein poly
    Vec<Vec<Bernstein<T> > > b(2,Vec<Bernstein<T> >(2,Bernstein<T>()));

    for (size_t i = 0; i<poly.size(); i++) {
        Vec<Bernstein<T> > bi = poly[i].splitAt(dim,split_point);
        for (size_t j = 0; j<bi.size(); j++)
            b[j][i] = bi[j];
    }

    Vec<Bernstein2d<T> > ret;
    for (size_t i = 0; i<2; i++) {
        ret.push_back(Bernstein2d<T>(b[i],rest));
    }
    return ret;
}

template <class T>
Vec<double> Bernstein2d<T>::dx() const
{
    Vec<double> dx = poly[0].dx();
    Vec<double> dy = poly[1].dx();
    return ((dx*dx)+(dy*dy));
}

template <class T>
size_t Bernstein2d<T>::maxDxDir() const
{
    if (isEmpty())
        throw("Bernstein2d::maxDxDir() called on empty Bernstein2d!");

    int max_idx;
    dx().max(&max_idx);

    if (max_idx < 0)
        throw("Bernstein2d::maxDxDir() would return idx < 0");

    return (size_t)max_idx;
}

template <class T>
void Bernstein2d<T>::setupPoints() const
{
    if (!points.empty()) {
        std::cout << " Bernstein2d::setupPoints - points not empty -> not computing!" << std::endl;
        return;
    }

    points = Polygon(poly[0].coeffs(),poly[1].coeffs());
}

template <>
void Bernstein2d<Iv>::setupPoints() const
{
    if (!points.empty()) {
        std::cout << " Bernstein2d::setupPoints - points not empty -> not computing!" << std::endl;
        return;
    }

    iv_points = IvPolygon(poly[0].coeffs(),poly[1].coeffs());
}

template <class T>
void Bernstein2d<T>::computeConvexHull() const
{
    if (!hull.empty()) {
        std::cout << " Bernstein2d::computeConvexHull - hull not empty -> not computing!" << std::endl;
        return;
    }

    if (points.empty() && iv_points.empty())
        setupPoints();

    if (!points.empty())
        hull = points.convexHull();
    else if (!iv_points.empty())
        hull = iv_points.convexHull();
}

template <class T>
void Bernstein2d<T>::storeConvexHull() const
{
    double t=sec();

    std::ofstream hnd("b_coeffs.m");
    if (!hnd.good()) {
        throw("Bernstein2d::storeConvexHull - Error opening output file!");
    }

    for (size_t i=0; i<poly.size(); i++) {

        hnd << "lims(" << (i+1) << ",:) = [ " << poly[i].bound().lower() << ", " << poly[i].bound().upper() << "];" << std::endl;
        hnd << "rest(" << (i+1) << ",:) = [ " << rest[i].lower() << ", " << rest[i].upper() << "];" << std::endl;
        hnd << "coeffs(" << (i+1) << ",:) = " << poly[i].coeffs() << ";" << std::endl;
    }

    hnd << "hull = " << hull << ";" << std::endl;

    hnd.close();

    double t2=sec();
    std::cout << "Store to file: " << t2-t << "sec!" << std::endl;
}

/*!
 * This algorithm evaluates zero-inclusion using the algorithm presented in [Schauss, Peer, Buss, TAC 2017] for Taylor
 * Models. It makes use of Bernstein2dPatch to check whether zero is included in a face of the parameter box.
 * This method is based on the algorithm introduced in  \cite Zettler1998.
 */
template <class T>
void Bernstein2d<T>::checkPatches() const
{
    Vec<uint8_t> orders = poly[0].getOrders();

    //two dimensions are swept resulting in 4 edges
    for(size_t e1=0; e1 < orders.size()-1; e1++) {
        if (orders[e1] == 0)
            continue;

        for(size_t e2=e1+1; e2 < orders.size(); e2++) {
            if (orders[e2] == 0)
                continue;

            //all but two dimensions are set to fixed values (0,order+1), for all possibilities
            // there are 2^(dim-2) possibilities
            for (size_t S_0=0; S_0<pow(2,orders.size()-2); S_0++) {
                ldiv_t divresult;
                divresult.rem = S_0;
                Vec<uint8_t> id(orders.size(),0);
                bool max_order_zero = false;
                for(size_t i=0; i<orders.size();i++) {
                    if ( (i==e1) || (i==e2) )
                        continue;

                    divresult = ldiv(divresult.rem, 2);
                    if (divresult.rem != 0) {
                        if (orders[i] == 0) {
                            max_order_zero = true;
                            break;
                        } else
                            id[i] = orders[i];
                    }

                    divresult.rem = divresult.quot;
                }
                if (max_order_zero)
                    continue;

                Vec<Polygon> P;

                Polygon pol;
                id[e2] = 0;
                for (int idx = 0; idx <= (int)orders[e1]; idx++) {
                    id[e1] = idx;
                    pol.push_back(points[poly[0].getIdx(id)]);
                }
                P.push_back(pol);

                pol.clear();
                id[e1] = orders[e1];
                for (int idx = 0; idx <= (int)orders[e2]; idx++) {
                    id[e2] = idx;
                    pol.push_back(points[poly[0].getIdx(id)]);
                }
                P.push_back(pol);

                pol.clear();
                id[e2] = orders[e2];
                for (int idx = (int)orders[e1]; idx >= 0; idx--) {
                    id[e1] = idx;
                    pol.push_back(points[poly[0].getIdx(id)]);
                }
                P.push_back(pol);

                pol.clear();
                id[e1] = 0;
                for (int idx = (int)orders[e2]; idx >= 0; idx--) {
                    id[e2] = idx;
                    pol.push_back(points[poly[0].getIdx(id)]);
                }
                P.push_back(pol);

                Bernstein2dPatch<Polygon> b2d_patch(P);
                if ( b2d_patch.includes(-rest[0],-rest[1]) ) {
                    if (b2d_patch.edgesInclude(-rest[0],-rest[1])) {
                        rest_partially_included = true;
                    } else {
                        rest_included = true;
                        rest_partially_included = false;
                        zero_included = false;  //TODO: should be true!
                        return;
                    }
                }

                if ( !rest_partially_included && b2d_patch.includesZero() ) {
                    if (!b2d_patch.edgesIncludeZero()) {
                        zero_included = true;
                    }
                }
            }
        }
    }
}

/*!
 * See checkPatches() for details, the only difference is that we consider Interval Bernstein coefficients in this
 * case.
 * \todo We should combine checkPatches() and checkPatchesIv(). To do this nicely we need a new class (simplePoly or
 * convexHull or so) is a template class for Iv or double and offers the functionality of IvPolygon.
 */
template <class T>
void Bernstein2d<T>::checkPatchesIv() const
{
    Vec<uint8_t> orders = poly[0].getOrders();

    //two dimensions are swept resulting in 4 edges
    for(size_t e1=0; e1 < orders.size()-1; e1++) {
        if (orders[e1] == 0)
            continue;

        for(size_t e2=e1+1; e2 < orders.size(); e2++) {
            if (orders[e2] == 0)
                continue;

            //all but two dimensions are set to fixed values (0,order+1), for all possibilities
            // there are 2^(dim-2) possibilities
            for (size_t S_0=0; S_0<pow(2,orders.size()-2); S_0++) {
                ldiv_t divresult;
                divresult.rem = S_0;
                Vec<uint8_t> id(orders.size(),0);
                bool max_order_zero = false;
                for(size_t i=0; i<orders.size();i++) {
                    if ( (i==e1) || (i==e2) )
                        continue;

                    divresult = ldiv(divresult.rem, 2);
                    if (divresult.rem != 0) {
                        if (orders[i] == 0) {
                            max_order_zero = true;
                            break;
                        } else
                            id[i] = orders[i];
                    }

                    divresult.rem = divresult.quot;
                }
                if (max_order_zero)
                    continue;

                Vec<IvPolygon> P;

                IvPolygon pol;
                id[e2] = 0;
                for (int idx = 0; idx <= (int)orders[e1]; idx++) {
                    id[e1] = idx;
                    pol.push_back(iv_points[poly[0].getIdx(id)]);
                }
                P.push_back(pol);

                pol.clear();
                id[e1] = orders[e1];
                for (int idx = 0; idx <= (int)orders[e2]; idx++) {
                    id[e2] = idx;
                    pol.push_back(iv_points[poly[0].getIdx(id)]);
                }
                P.push_back(pol);

                pol.clear();
                id[e2] = orders[e2];
                for (int idx = (int)orders[e1]; idx >= 0; idx--) {
                    id[e1] = idx;
                    pol.push_back(iv_points[poly[0].getIdx(id)]);
                }
                P.push_back(pol);

                pol.clear();
                id[e1] = 0;
                for (int idx = (int)orders[e2]; idx >= 0; idx--) {
                    id[e2] = idx;
                    pol.push_back(iv_points[poly[0].getIdx(id)]);
                }
                P.push_back(pol);

                Bernstein2dPatch<IvPolygon> b2d_patch(P);
                if ( b2d_patch.includes(-rest[0],-rest[1]) ) {
                    if (b2d_patch.edgesInclude(-rest[0],-rest[1])) {
                        rest_partially_included = true;
                    } else {
                        rest_included = true;
                        rest_partially_included = false;
                        zero_included = false;  //TODO: should be true!
                        return;
                    }
                }

                if ( !rest_partially_included && b2d_patch.includesZero() ) {
                    if (!b2d_patch.edgesIncludeZero()) {
                        zero_included = true;
                    }
                }
            }
        }
    }
}

template class Bernstein2d<double>;
template class Bernstein2d<Iv>;

template <class T>
std::ostream & operator<<(std::ostream &os, const Bernstein2d<T>& b2d)
{
    os << "hull:          " << b2d.hull;
    os << "rest:          " << b2d.rest << std::endl;
    os << "restIncluded:          " << b2d.restIncluded() << std::endl;
    os << "restExcluded:          " << b2d.restExcluded() << std::endl;
    os << "restPartiallyIncluded: " << b2d.restPartiallyIncluded() << std::endl;
    os << "restPartiallyExcluded: " << b2d.restPartiallyExcluded() << std::endl;
    os << "zeroIncluded:          " << b2d.zeroIncluded() << std::endl;
    os << "zeroExcluded:          " << b2d.zeroExcluded() << std::endl;
    os << "num_points:            " << b2d.points.size() << std::endl;
    return os;
}

template std::ostream & operator<<(std::ostream &os, const Bernstein2d<double>& b2d);
template std::ostream & operator<<(std::ostream &os, const Bernstein2d<Iv>& b2d);
