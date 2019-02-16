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

#include "bernstein_2d_patch.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "helpers.h"
#include "iv_polygon.h"

template <class T>
unsigned int Bernstein2dPatch<T>::count = 0;

template <class T>
Bernstein2dPatch<T>::Bernstein2dPatch(Vec<T> &edges)
{
    id = count;
    count++;

    for (size_t i=0; i<edges.size(); i++) {
        R.push_back(edges[i][0]);
    }
    R.touchesIncludes(true);
    B = edges;
}

template <>
Bernstein2dPatch<IvPolygon>::Bernstein2dPatch(Vec<IvPolygon> &edges)
{
    id = count;
    count++;

    for (size_t i=0; i<edges.size(); i++) {
        R.push_back(Polygon::Point_2(edges[i][0][0].mid(),edges[i][0][1].mid()));
    }
    R.touchesIncludes(true);
    B = edges;
}

/*!
 * \todo Check whether we can remove R.includesZero as this is not self explanatory (we could, e.g., only do that fast
 * check if zero is included in the Intervals)
 * \todo Rename fuction as we would expect this to check inclusion in the face and not in the patch connecting the
 * corners (actually we should keep an include function which checks R.includes && !edgesInclude)
 */
template <class T>
bool Bernstein2dPatch<T>::includes(const Iv &x, const Iv &y) const{
    return ( (R.includesZero()) && (R.includes(x,y) == Polygon::COMPLETELY_INCLUDED));
}

/*!
 * \todo Rename fuction as we would expect this to check inclusion in the face and not in the patch connecting the
 * corners (actually we should keep an include function which checks R.includes && !edgesInclude)
 */
template <class T>
bool Bernstein2dPatch<T>::includesZero() const{
    return R.includesZero();
}

template <class T>
bool Bernstein2dPatch<T>::edgesInclude(const Iv &x, const Iv &y) const
{
    for (size_t i = 0; i<B.size(); i++) {
        Polygon C(B[i].convexHull());
        //C.touchesExcludes(true);
        Polygon::include_t b = C.includes(x,y);
        if ( (b==Polygon::COMPLETELY_INCLUDED) || (b==Polygon::PARTIALLY_INCLUDED) ) {
            return true;
            break;
        }
    }
    return false;
}

template <class T>
bool Bernstein2dPatch<T>::edgesIncludeZero() const
{
    for (size_t i = 0; i<B.size(); i++) {
        Polygon C(B[i].convexHull());
        //C.touchesExcludes(true);
        if ( C.includesZero() ) {
            return true;
            break;
        }
    }
    return false;
}

template <class T>
Bernstein2dPatch<T>::~Bernstein2dPatch()
{
}

template class Bernstein2dPatch<Polygon>;
template class Bernstein2dPatch<IvPolygon>;

template <class T>
std::ostream & operator<<(std::ostream &os, const Bernstein2dPatch<T>& b)
{
    for (size_t i=0; i<b.B.size(); i++) {
        os<< "edges(" << (b.id+1) << ").B" << i << " = " << b.B[i] << std::endl;
    }
    os << "edges(" << (b.id+1) << ").R = " << b.R << std::endl;

    return os;
}

template std::ostream & operator<<(std::ostream &os, const Bernstein2dPatch<Polygon>& b);
template std::ostream & operator<<(std::ostream &os, const Bernstein2dPatch<IvPolygon>& b);
