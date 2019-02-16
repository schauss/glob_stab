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

#define BOOST_TEST_MODULE polynomial

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "polynomial.h"
#include "bernstein.h"
#include <iostream>

BOOST_AUTO_TEST_CASE(simple_polynomial)
{
    // create empty polynomial
    Polynomial<double> p;
    //BOOST_CHECK(p.empty());
    std::cout << "Empty polynomial: " << p << std::endl;

    Polynomial<double> a("a");
    std::cout << "a:" << a << std::endl;
    Polynomial<double> b("b");
    std::cout << "b:" << b << std::endl;

    Polynomial<double> c("c");
    std::cout << "c:" << c << std::endl;

    std::cout << (-7.654+3*a*b-0.5*a*a)*(c*c+2) << std::endl;
}

BOOST_AUTO_TEST_CASE(fixed_polynomial)
{
    // Fix variables
    Vec<std::string> vars;
    vars.push_back("a");
    vars.push_back("b");
    vars.push_back("c");
    Monomial<double>::fixVariables(vars);

    // create empty polynomial
    Polynomial<double> p;
    //BOOST_CHECK(p.empty());
    std::cout << "Empty polynomial: " << p << std::endl;

    Polynomial<double> a("a");
    std::cout << "a:" << a << std::endl;
    Polynomial<double> b("b");
    std::cout << "b:" << b << std::endl;

    Polynomial<double> c("c");
    std::cout << "c:" << c << std::endl;

    std::cout << (-7.654+3*a*b-0.5*a*a)*(c*c+2) << std::endl;
}

BOOST_AUTO_TEST_CASE(shift_polynomial)
{
    // Fix variables
    Vec<std::string> vars;
    vars.push_back("a");
    vars.push_back("b");
    vars.push_back("c");
    Monomial<double>::fixVariables(vars);

    Polynomial<double> a("a");
    Polynomial<double> b("b");
    Polynomial<double> c("c");

    Polynomial<double> p = (2*a -3*b+4*c)*(-3*a*b+7*c*c*c+2*a*a*b-17*b*b*c);
    std::cout << "p = " << p << std::endl;

    Polynomial<double> q = p.changeDomain(Iv(0,1),Iv(-1,1));
    std::cout << "q = p shifted from [0,1] to [-1,1] = " << q << std::endl;

    std::cout << "p.bound(Iv(0,1)) = " << p.bound(Iv(0,1)) << std::endl;
    std::cout << "q.bound(Iv(-1,1)) = " << q.bound(Iv(-1,1)) << std::endl;

    Polynomial<double> r = q.changeDomain(Iv(-1,1),Iv(0,1));
    std::cout << "r = q shifted from [-1,1] to [0,1] = p = " << r << std::endl;
    std::cout << "r.bound(Iv(0,1)) = " << r.bound(Iv(0,1)) << std::endl;
}

BOOST_AUTO_TEST_CASE(bernstein_transformation)
{
    // Fix variables
    Vec<std::string> vars;
    vars.push_back("a");
    vars.push_back("b");
    vars.push_back("c");
    Monomial<double>::fixVariables(vars);

    Polynomial<double> a("a");
    Polynomial<double> b("b");
    Polynomial<double> c("c");

    Polynomial<double> p = (2*a -3*b+4*c)*(-3*a*b+7*c*c*c+2*a*a*b-17*b*b*c);
    std::cout << std::endl << "p = " << p << std::endl;

    std::cout << std::endl << "Bernstein-bound: " << p.bound(Iv(0,1),true) << std::endl;
    std::cout << std::endl << "Iv-bound: " << p.bound(Iv(0,1),false) << std::endl;

    p = (6*b+3*a+4*c)*c;
    std::cout << std::endl << "p = " << p << std::endl;

    Bernstein<double> B(p);
    std::cout << std::endl << "Bernstein polynomial B: " << std::endl << B << std::endl;
}

BOOST_AUTO_TEST_CASE(simple_polynomial_Iv)
{
    // create empty polynomial
    Polynomial<Iv> p;
    //BOOST_CHECK(p.empty());
    std::cout << "Empty polynomial: " << p << std::endl;

    Polynomial<Iv> a("a");
    std::cout << "a:" << a << std::endl;
    Polynomial<Iv> b("b");
    std::cout << "b:" << b << std::endl;

    Polynomial<Iv> c("c");
    std::cout << "c:" << c << std::endl;

    std::cout << (-7.654+3*a*b-0.5*a*a)*(c*c+2) << std::endl;
}

BOOST_AUTO_TEST_CASE(fixed_polynomial_Iv)
{
    // Fix variables
    Vec<std::string> vars;
    vars.push_back("a");
    vars.push_back("b");
    vars.push_back("c");
    Monomial<Iv>::fixVariables(vars);

    // create empty polynomial
    Polynomial<Iv> p;
    //BOOST_CHECK(p.empty());
    std::cout << "Empty polynomial: " << p << std::endl;

    Polynomial<Iv> a("a");
    std::cout << "a:" << a << std::endl;
    Polynomial<Iv> b("b");
    std::cout << "b:" << b << std::endl;

    Polynomial<Iv> c("c");
    std::cout << "c:" << c << std::endl;

    std::cout << (-7.654+3*a*b-0.5*a*a)*(c*c+2) << std::endl;
}

BOOST_AUTO_TEST_CASE(shift_polynomial_Iv)
{
    // Fix variables
    Vec<std::string> vars;
    vars.push_back("a");
    vars.push_back("b");
    vars.push_back("c");
    Monomial<Iv>::fixVariables(vars);

    Polynomial<Iv> a("a");
    Polynomial<Iv> b("b");
    Polynomial<Iv> c("c");

    Polynomial<Iv> p = (2*a -3*b+4*c)*(-3*a*b+7*c*c*c+2*a*a*b-17*b*b*c);
    std::cout << "p = " << p << std::endl;

    Polynomial<Iv> q = p.changeDomain(Iv(0,1),Iv(-1,1));
    std::cout << "q = p shifted from [0,1] to [-1,1] = " << q << std::endl;

    std::cout << "p.bound(Iv(0,1)) = " << p.bound(Iv(0,1)) << std::endl;
    std::cout << "q.bound(Iv(-1,1)) = " << q.bound(Iv(-1,1)) << std::endl;

    Polynomial<Iv> r = q.changeDomain(Iv(-1,1),Iv(0,1));
    std::cout << "r = q shifted from [-1,1] to [0,1] = p = " << r << std::endl;
    std::cout << "r.bound(Iv(0,1)) = " << r.bound(Iv(0,1)) << std::endl;
}

BOOST_AUTO_TEST_CASE(bernstein_transformation_Iv)
{
    // Fix variables
    Vec<std::string> vars;
    vars.push_back("a");
    vars.push_back("b");
    vars.push_back("c");
    Monomial<Iv>::fixVariables(vars);

    Polynomial<Iv> a("a");
    Polynomial<Iv> b("b");
    Polynomial<Iv> c("c");

    Polynomial<Iv> p = (2*a -3*b+4*c)*(-3*a*b+7*c*c*c+2*a*a*b-17*b*b*c);
    std::cout << std::endl << "p = " << p << std::endl;

    std::cout << std::endl << "Bernstein-bound: " << p.bound(Iv(0,1),true) << std::endl;
    std::cout << std::endl << "Iv-bound: " << p.bound(Iv(0,1),false) << std::endl;

    p = (6*b+3*a+4*c)*c;
    std::cout << std::endl << "p = " << p << std::endl;

    Bernstein<Iv> B(p);
    std::cout << std::endl << "Bernstein polynomial B: " << std::endl << B << std::endl;
}
