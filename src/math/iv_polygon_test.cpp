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

#define BOOST_TEST_MODULE iv_polygon

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "iv.h"
#include <iostream>

void line_intersection(Iv x1, Iv y1, Iv x2, Iv y2, Iv x3, Iv y3, Iv x4, Iv y4)
{
    Iv x12 = x1-x2;
    Iv x34 = x3-x4;
    Iv y12 = y1-y2;
    Iv y34 = y3-y4;

    Iv den = x12*y34 - y12*x34;

    //std::cout << "den = " << den << std::endl;

    Iv xy12 = x1*y2-y1*x2;
    Iv xy34 = x3*y4-y3*x4;

    Iv x = (xy12*x34-xy34*x12)/den;
    Iv y = (xy12*y34-xy34*y12)/den;

    //std::cout << "x = " << x << std::endl;
    //std::cout << "y = " << y << std::endl;

    if ( ( (x.lower() > x1.upper()) && (x.lower() > x2.upper()) ) ||
         ( (x.upper() < x1.lower()) && (x.upper() < x2.lower()) ) )
        std::cout << "-> No intersection!" << std::endl;
    else if ( ( (x.lower() > x1.upper()) && (x.upper() < x2.lower()) ) ||
         ( (x.lower() > x2.upper()) && (x.upper() < x1.lower()) ) )
        std::cout << "-> Cetain intersection!" << std::endl;
    else
        std::cout << "-> Unknown!" << std::endl;
}

void line_intersection_1(Iv x1, Iv y1, Iv x2, Iv y2, Iv x3, Iv y3, Iv x4, Iv y4)
{
    Vec<Iv> in;
    in.push_back(x1);
    in.push_back(y1);
    in.push_back(x2);
    in.push_back(y2);
    in.push_back(x3);
    in.push_back(y3);
    in.push_back(x4);
    in.push_back(y4);

    Vec<double> out(8,0.0);

    for (size_t i=0; i<256; i++) {
        ldiv_t divresult;
        divresult.rem = i;
        for(size_t j=0; j<8;j++) {
            divresult = ldiv(divresult.rem, 2);
            if (divresult.rem != 0) {
                std::cout << "1";
                out[j] = in[j].upper();
            } else {
                std::cout << "0";
                out[j] = in[j].lower();
            }

            divresult.rem = divresult.quot;
        }
        std::cout << ": ";
        line_intersection(out[0],out[1],out[2],out[3],out[4],out[5],out[6],out[7]);
    }
}

BOOST_AUTO_TEST_CASE(iv_polygon)
{
    Iv x1(0,2);
    Iv y1(0.0);

    Iv x2(1.0);
    Iv y2(1.0);

    Iv x3(0.0);
    Iv y3(1.0);

    Iv x4(1.0);
    Iv y4(0.0);

    line_intersection(x1,y1,x2,y2,x3,y3,x4,y4);
    line_intersection_1(x1,y1,x2,y2,x3,y3,x4,y4);
}
