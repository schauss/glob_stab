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

#define BOOST_TEST_MODULE TaylorModel

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "taylor_model.h"
#include <iostream>
#include <iomanip>

BOOST_AUTO_TEST_CASE(simple_TaylorModel)
{
    // create empty polynomial
    TaylorModel<> p;
    //BOOST_CHECK(p.empty());

    TaylorModel<> a("a");
    TaylorModel<> b("b");
    TaylorModel<> c("c");

    std::cout << (-7.654+3*a*b-0.5*a*a)*(c*c+2) << std::endl;
}

BOOST_AUTO_TEST_CASE(fixed_TaylorModel)
{
    // Fix variables
    Vec<std::string> vars;
    vars.push_back("a");
    vars.push_back("b");
    vars.push_back("c");
    Monomial<double>::fixVariables(vars);

    // create empty polynomial
    TaylorModel<> p;
    //BOOST_CHECK(p.empty());

    TaylorModel<> a("a");
    TaylorModel<> b("b");
    TaylorModel<> c("c");

    std::cout << std::scientific << std::setprecision(15) << (-7.654+3*a*b-0.5*a*a)*(c*c+2) << std::endl;

    std::cout << std::scientific << std::setprecision(15) << "div: " << (2.0 + a*2.0)/(1.0+b*0.001)*(-7.654+3*a*b-0.5*a*a)*(c*c+2)*(-7.654+3*a*b-0.5*a*a)*(c*c+2) << std::endl;
}

BOOST_AUTO_TEST_CASE(simple_TaylorModel_Iv)
{
    // create empty polynomial
    TaylorModel<Iv> p;
    //BOOST_CHECK(p.empty());

    TaylorModel<Iv> a("a");
    TaylorModel<Iv> b("b");
    TaylorModel<Iv> c("c");

    std::cout << (-7.654+3*a*b-0.5*a*a)*(c*c+2) << std::endl;
}

BOOST_AUTO_TEST_CASE(fixed_TaylorModel_Iv)
{
    // Fix variables
    Vec<std::string> vars;
    vars.push_back("a");
    vars.push_back("b");
    vars.push_back("c");
    Monomial<Iv>::fixVariables(vars);

    // create empty polynomial
    TaylorModel<Iv> p;
    //BOOST_CHECK(p.empty());

    TaylorModel<Iv> a("a");
    TaylorModel<Iv> b("b");
    TaylorModel<Iv> c("c");

    std::cout << std::scientific << std::setprecision(15) << (-7.654+3*a*b-0.5*a*a)*(c*c+2) << std::endl;

    std::cout << std::scientific << std::setprecision(15) << "div: " << (2.0 + a*2.0)/(1.0+b*0.001)*(-7.654+3*a*b-0.5*a*a)*(c*c+2)*(-7.654+3*a*b-0.5*a*a)*(c*c+2) << std::endl;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "sin(a) = " << sin(a) << std::endl;
    std::cout << "sin(10*a) = " << sin(10*a) << std::endl;

    std::cout << "exp(a) = " << exp(a) << std::endl;
    std::cout << "exp(10*a) = " << exp(10*a) << std::endl;
    std::cout << "exp(9.9+0.1*a) = " << exp(9.9+0.1*a) << std::endl;
}

BOOST_AUTO_TEST_CASE(special)
{
    // Fix variables
    Vec<std::string> vars;
    vars.push_back("sigma");
    vars.push_back("omega");
    vars.push_back("td");
    Monomial<Iv>::fixVariables(vars);

    std::cout << std::endl << "##########################################" << std::endl;
    std::cout << " BOARDER " << std::endl;
    std::cout << "##########################################" << std::endl << std::endl;

    //TaylorModel<Iv> sigma(0.0+TaylorModel<Iv>("sigma"));
    TaylorModel<Iv> sigma(0.0);
    TaylorModel<Iv> omega(0.0+1*TaylorModel<Iv>("omega"));

    std::complex<TaylorModel<Iv> > s(false,omega);
    std::cout << std::endl << "s = " << s << std::endl;

    std::complex<TaylorModel<Iv> > td(0.1,false);
    std::cout << std::endl << "td = " << td << std::endl;

    std::cout << std::endl << "exp(-td*s) = " << exp(-td*s) << std::endl;

    std::complex<TaylorModel<Iv> > s_(false,omega.inv());
    std::cout << std::endl << "s_ = " << s_ << std::endl;

    std::cout << std::endl << "(-td*s_) = " << (-td*s_) << std::endl;

    std::cout << std::endl << "exp(-td*s_) = " << exp(-td*s_) << std::endl;

    std::cout << std::endl << "##########################################" << std::endl;
    std::cout << " RIGHT HALF PLANE " << std::endl;
    std::cout << "##########################################" << std::endl << std::endl;

    sigma = TaylorModel<Iv>(0.0+TaylorModel<Iv>("sigma"));

    std::complex<TaylorModel<Iv> > s1(sigma,omega);
    std::cout << std::endl << "s1 = " << s1 << std::endl;

    std::cout << std::endl << "td = " << td << std::endl;

    std::cout << std::endl << "exp(-td*s1) = " << exp(-td*s1) << std::endl;

    std::complex<TaylorModel<Iv> > s2(sigma,omega.inv());
    std::cout << std::endl << "s2 = " << s2 << std::endl;

    std::cout << std::endl << "(-td*s2) = " << (-td*s2) << std::endl;

    std::cout << std::endl << "exp(-td*s2) = " << exp(-td*s2) << std::endl;

    std::complex<TaylorModel<Iv> > s3(sigma.inv(),omega);
    std::cout << std::endl << "s3 = " << s3 << std::endl;

    std::cout << std::endl << "(-td*s3) = " << (-td*s3) << std::endl;

    std::cout << std::endl << "exp(-td*s3) = " << exp(-td*s3) << std::endl;

    std::complex<TaylorModel<Iv> > s4(sigma.inv(),omega.inv());
    std::cout << std::endl << "s4 = " << s4 << std::endl;

    std::cout << std::endl << "(-td*s4) = " << (-td*s4) << std::endl;

    std::cout << std::endl << "exp(-td*s4) = " << exp(-td*s4) << std::endl;

    std::cout << std::endl << "##########################################" << std::endl;
    std::cout << " RIGHT HALF PLANE INVERSE (as implemented)" << std::endl;
    std::cout << "##########################################" << std::endl << std::endl;

    omega = TaylorModel<Iv>(0.0+1.0001*TaylorModel<Iv>("omega"));
    //omega = TaylorModel<Iv>(0.0001+0.1*TaylorModel<Iv>("omega"));
    //sigma = TaylorModel<Iv>(0.0+1.0001*TaylorModel<Iv>("sigma"));
    //sigma = TaylorModel<Iv>(0.50005+0.50005*TaylorModel<Iv>("sigma"));
    sigma = TaylorModel<Iv>(0.5+0.1*TaylorModel<Iv>("sigma"));

    s = std::complex<TaylorModel<Iv> >(sigma,omega);
    std::cout << std::endl << "s = " << s << std::endl;

    TaylorModel<Iv>::correctTrigonometric(true);

    std::complex<TaylorModel<Iv> > s_inv = std::complex<TaylorModel<Iv> >(1.0)/s;
    std::cout << std::endl << "s_inv = " << s_inv << std::endl;

    std::cout << std::endl << "(-td*s_inv) = " << (-td*s_inv) << std::endl;

    std::cout << std::endl << "exp(-td*s_inv) = " << exp(-td*s_inv) << std::endl;

//    std::cout << std::endl << "##########################################" << std::endl;
//    std::cout << " SIMPLER TESTS" << std::endl;
//    std::cout << "##########################################" << std::endl << std::endl;

//    std::cout << "TM(1)/sigma = " << TaylorModel<Iv>(1.0)/sigma << std::endl << std::endl;
//    std::cout << "1.0/sigma = " << 1.0/sigma << std::endl << std::endl;
//    std::cout << "(1.0/sigma)*sigma = " << (1.0/sigma)*sigma << std::endl << std::endl;

}
