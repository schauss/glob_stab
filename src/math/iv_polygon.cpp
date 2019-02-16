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

#include "iv_polygon.h"

IvPolygon::IvPolygon()
{

}

IvPolygon::IvPolygon(const Vec<Iv> &x, const Vec<Iv> &y)
    : x(x), y(y)
{
    if (x.size() != y.size())
        throw("IvPolygon::IvPolygon(x,y) - different sizes for x and y not allowed!");

}

IvPolygon::~IvPolygon()
{

}

/*!
 * The convex hull is determined by first creating a Polygon consisting of the
 * four corners of each interval-box and then calling Polygon::convexHull().
 */
Polygon IvPolygon::convexHull() const
{
    Vec<Polygon::Point_2> p;
    p.reserve(x.size()*4);
    for (size_t i=0; i<x.size();i++) {
        p.push_back(Polygon::Point_2(x[i].lower(),y[i].lower()));
        p.push_back(Polygon::Point_2(x[i].upper(),y[i].lower()));
        p.push_back(Polygon::Point_2(x[i].lower(),y[i].upper()));
        p.push_back(Polygon::Point_2(x[i].upper(),y[i].upper()));
    }

    return Polygon(p).convexHull();
}

/*!
 * The order of the points in the Polygon is the same as the order of the
 * interval-boxes in IvPolygon::x and IvPolygon::y.
 */
Polygon IvPolygon::mid() const
{
    Vec<Polygon::Point_2> p;
    p.reserve(x.size());
    for (size_t i=0; i<x.size();i++) {
        p.push_back(Polygon::Point_2(x[i].mid(),y[i].mid()));
    }

    return Polygon(p);
}

void IvPolygon::push_back(const Vec<Iv> p)
{
    if (p.size() != 2)
        throw("IvPolygon::push_back must be called with Vec of size 2!");

    x.push_back(p[0]);
    y.push_back(p[1]);
}


void IvPolygon::clear()
{
    x.clear();
    y.clear();
}

Vec<Iv> IvPolygon::operator[](size_t i) const
{
    if (i>=size())
        throw("IvPolygon::[] - bound violation!");

    Vec<Iv> ret;
    ret.push_back(x[i]);
    ret.push_back(y[i]);
    return ret;
}

std::ostream & operator<<(std::ostream &os, const IvPolygon& p)
{
    os << "[ ";
    for (size_t i=0; i<p.x.size(); i++) {
        os << "[ " << p.x[i] << "; " << p.y[i] << " ]" << (i==(p.x.size()-1)?" ]":", ");
    }

    return os;
}
