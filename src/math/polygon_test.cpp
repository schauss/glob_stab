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

#define BOOST_TEST_MODULE polygon

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "polygon.h"
#include <iostream>

BOOST_AUTO_TEST_CASE(simple_polygon)
{
    // create empty polygon
    Polygon p;
    BOOST_CHECK(p.empty());

    // add points and validate, that empty is working!
    p.push_back(Polygon::Point_2(0.0,0.0));
    BOOST_CHECK(p.empty());
    p.push_back(Polygon::Point_2(1.0,0.0));
    BOOST_CHECK(p.empty());
    p.push_back(Polygon::Point_2(1.0,1.0));
    BOOST_CHECK(!p.empty());
    p.push_back(Polygon::Point_2(0.0,1.0));
    BOOST_CHECK(!p.empty());

    // check inclusion of (double,double) for simple polygon
    BOOST_CHECK(p.includes(0.0,0.0)==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(1.0,0.3)==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(0.7,0.3)==Polygon::COMPLETELY_INCLUDED);
    BOOST_CHECK(p.includes(1.1,0.3)==Polygon::NOT_INCLUDED);

    // add fifth vertex to simple polygon
    p.push_back(Polygon::Point_2(0.0,0.5));
    BOOST_CHECK(!p.empty());
    BOOST_CHECK(p.includes(0.0,0.0)==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(1.0,0.3)==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(0.7,0.3)==Polygon::COMPLETELY_INCLUDED);
    BOOST_CHECK(p.includes(1.1,0.3)==Polygon::NOT_INCLUDED);

    //check interval inclusion
    BOOST_CHECK(p.includes(Iv(0.4,0.6),Iv(0.4,0.6))==Polygon::COMPLETELY_INCLUDED);
    BOOST_CHECK(p.includes(Iv(0.4,1.1),Iv(0.4,0.6))==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(Iv(0.9999,1.1),Iv(0.4,0.6))==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(Iv(1.01,1.1),Iv(0.4,0.6))==Polygon::NOT_INCLUDED);
    BOOST_CHECK(p.includes(Iv(1.0,1.1),Iv(0.4,0.6))==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(Iv(1.1,1.2),Iv(0.4,0.6))==Polygon::NOT_INCLUDED);
}

BOOST_AUTO_TEST_CASE(complex_polygon)
    {
    // create empty polygon
    Polygon p;
    BOOST_CHECK(p.empty());

    // make stupid collinear polygon and check empty-function -> not working, throws!
    /*
    p.clear();
    p.push_back(Polygon::Point_2(0.0,0.0));
    p.push_back(Polygon::Point_2(2.0,0.0));
    p.push_back(Polygon::Point_2(1.0,0.0));
    p.push_back(Polygon::Point_2(3.0,0.0));
    BOOST_CHECK(p.empty());
    */

    // make non-simple polygon and check empty-function
    p.clear();
    p.push_back(Polygon::Point_2(0.0,0.0));
    p.push_back(Polygon::Point_2(1.0,2.0));
    p.push_back(Polygon::Point_2(0.0,2.0));
    p.push_back(Polygon::Point_2(1.0,0.0));
    BOOST_CHECK(!p.empty());

    // check inclusion in non-simple polygon
    BOOST_CHECK(p.includes(0.5,0.0)==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(0.5,1.0)==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(0.5,1.1)==Polygon::COMPLETELY_INCLUDED);
    BOOST_CHECK(p.includes(0.5,0.1)==Polygon::COMPLETELY_INCLUDED);
    BOOST_CHECK(p.includes(0.5,2.1)==Polygon::NOT_INCLUDED);
    BOOST_CHECK(p.includes(0.6,1.0)==Polygon::NOT_INCLUDED);

    //check interval inclusion
    BOOST_CHECK(p.includes(Iv(0.4,0.6),Iv(0.01,0.2))==Polygon::COMPLETELY_INCLUDED);
    BOOST_CHECK(p.includes(Iv(0.4,0.6),Iv(0.0,0.2))==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(Iv(0.4,0.6),Iv(-0.2,0.2))==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(Iv(0.4,0.6),Iv(0.9,1.1))==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(Iv(0.4,0.6),Iv(2.01,2.1))==Polygon::NOT_INCLUDED);
    BOOST_CHECK(p.includes(Iv(0.4,0.6),Iv(2.0,2.1))==Polygon::PARTIALLY_INCLUDED);
    BOOST_CHECK(p.includes(Iv(0.4,0.6),Iv(2.1,2.2))==Polygon::NOT_INCLUDED);

    // add fifth vertex to non-simple polygon -> not working, throws!
    /*
    p.push_back(Polygon::Point_2(0.5,0.0));
    BOOST_CHECK(!p.empty());
    BOOST_CHECK(p.includes(0.5,1.0)==Polygon::PARTIALLY_INCLUDED);
    */
}
