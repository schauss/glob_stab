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

#define BOOST_TEST_MODULE ndArray

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "nd_array.hpp"
#include <iostream>

BOOST_AUTO_TEST_CASE(simple_tests)
{
    std::vector<size_t> dims;
    dims.push_back(3);
    dims.push_back(4);
    dims.push_back(5);

    ndArray<double> a(dims);
    a.print();

    std::cout << "Setting ndArray[1][2] = 5\n";
    a[1][2] = 5;

    std::vector<size_t> start, end;
    start.push_back(2);
    end.push_back(3);
    start.push_back(2);
    end.push_back(4);
    start.push_back(3);
    end.push_back(5);

    std::cout << "Setting ndArray([2,2,3],[3,4,5]) = 3\n";
    a(start,end) = 3;

    BOOST_CHECK(a(start,end)==3.0);
    BOOST_CHECK(!(a(start,end)==4.0));

    BOOST_CHECK(a[1][2]==5.0);
    BOOST_CHECK(a[0]==0.0);
    BOOST_CHECK(a[1][0]==0.0);
    BOOST_CHECK(!(a[1]==0.0));
    BOOST_CHECK(!(a[2]==0.0));

    start[0]=1;
    BOOST_CHECK(!(a(start,end)==3.0));

    for(size_t i=0; i<dims[0]; i++) {
        for(size_t j=0; j<dims[1]; j++) {
            for(size_t k=0; k<dims[2]; k++) {
                //std::cout << "a[" << i << "][" << j << "][" << k << "] = " << a[i][j][k] << "\n";

/*
                if ( (i==1) && (j==2) )
                    BOOST_CHECK( a[i][j][k] == 5.0 );
                else
                    BOOST_CHECK( a[i][j][k] == 0.0 );
*/
            }
        }
    }
}
