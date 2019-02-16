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

#define BOOST_TEST_MODULE iv

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "iv_sse2.h"
#include "iv.h"
#include <iostream>

BOOST_AUTO_TEST_CASE(iv_sse2)
{
    // create empty polygon
    IvSSE2<> a(1.23456789012345,1.23456789012345);
    std::cout << "a = " << a << std::endl;

    IvSSE2<> b(1,2);

    for (size_t i=1; i<1000000; i++) {
        b += a;
    }

    for (size_t i=1; i<1000000; i++) {
        b -= a;
    }

    BOOST_CHECK( (b.lower() < 1) && (b.upper() > 2));
    std::cout << "b = " << b << std::endl;

    b = IvSSE2<>(1,2);

    b += 1e5;
    for (size_t i=1; i<1000000; i++) {
        b += a;
        b -= a;
    }
    b -= 1e5;

    BOOST_CHECK( (b.lower() < 1) && (b.upper() > 2));
    std::cout << "b = " << b << std::endl;

    b = a;

    for (size_t i=1; i<1000000; i++) {
        b *= a;
        b /= a;
    }

    b -= a;

    BOOST_CHECK( (b.lower() < 0) && (b.upper() > 0));

    std::cout << "b = " << b << std::endl;
    //BOOST_CHECK(p.empty());
}

BOOST_AUTO_TEST_CASE(iv_trigonometric)
{
    Iv a;

    double offset = -5;
    double step   = 0.1;

    std::cout << std::endl;

    for (int i=0; i<100; i++) {
        a = Iv(offset,offset+step*i);
        std::cout << "sin(" << a << ") = " << a.sin() << std::endl;
    }

    std::cout << std::endl;

    for (int i=0; i<100; i++) {
        a = Iv(offset,offset+step*i);
        std::cout << "cos(" << a << ") = " << a.cos() << std::endl;
    }

    std::cout << std::endl;

    for (int i=0; i<100; i++) {
        a = Iv(offset,offset+step*i);
        std::cout << "exp(" << a << ") = " << a.exp() << std::endl;
    }
}
