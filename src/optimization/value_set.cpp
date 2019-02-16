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

#include "value_set.h"

#include "equation.h"
#include "parameters.h"
#include "taylor_model.h"
#include "config.h"

#include <boost/math/special_functions/binomial.hpp>

ValueSet::ValueSet(const Equation &chEq, double delta)
    : ch_eq(chEq),
      boundary(false),
      complex(false),
      inverse(false),
      delta(delta)
{

}

template <class T>
Vec<T> ValueSet::value(const std::map<std::string,T> &param) const
{
    T omega = get(param,"omega");
    T sigma;
    if (boundary)
        sigma = 0.0;
    else
        sigma = get(param,"sigma");

    std::ostringstream buf;

    if (CONFIG_VAR(verbosity,int) > 2) {
        buf << "##########################################" << std::endl;
        buf << "Evaluating value set for:" << std::endl << std::endl;
        buf << "boundary = " << boundary << std::endl;
        buf << "inverse = " << inverse << std::endl;
        buf << "complex = " << complex << std::endl;
        buf << "delta = " << delta << std::endl << std::endl;
        buf << "omega = " << omega << std::endl;
        buf << "sigma = " << sigma << std::endl << std::endl;
    }

    Vec<std::complex<T> > chEq;

    if (complex) {
        std::map<std::string,std::complex<T> > p;
        for (typename std::map<std::string,T>::const_iterator it=param.begin(); it!=param.end(); it++)
            p[it->first] = std::complex<T>(it->second,false);
        p["s"] = s(sigma,omega);
        if (CONFIG_VAR(verbosity,int) > 2)
            buf << "s = " << p["s"] << std::endl << std::endl;
        chEq = ch_eq.value(p);
    } else
        chEq = ch_eq.value(param);

    if (CONFIG_VAR(verbosity,int) > 2)
        buf << "chEq = " << chEq << std::endl << std::endl;

    if (delta != 0.0)
        sigma -= delta*omega;

    std::complex<T> eq(false,false);

    for (size_t i=0; i<chEq.size(); i++) {
        std::complex<T> eq_tmp(false,false);
        eq_tmp = s(sigma,omega,i,chEq.size()-1)*chEq[i];
        if (CONFIG_VAR(verbosity,int) > 2) {
            buf << "chEq[" << i << "] = " << eq_tmp << std::endl << std::endl;
        }
        eq += eq_tmp;
    }

    if (CONFIG_VAR(verbosity,int) > 2) {
        buf << "vs = " << eq << std::endl << std::endl;
        buf << "##########################################" << std::endl;
        std::cout << buf.str() << std::endl;
    }

    Vec<T> res(2,0.0);
    res[0] = eq.real();
    res[1] = eq.imag();

    return res;
}

template Vec<TaylorModel<double> > ValueSet::value<TaylorModel<double> >(const std::map<std::string,TaylorModel<double> > &param) const;
template Vec<TaylorModel<Iv> > ValueSet::value<TaylorModel<Iv> >(const std::map<std::string,TaylorModel<Iv> > &param) const;

template <class T>
T ValueSet::valueS0(const std::map<std::string,T> &param) const
{
    Equation ch_eq0 = ch_eq;
    ch_eq0.set("s",0.0);

    return ch_eq0.value(param)[0];
}

template TaylorModel<double> ValueSet::valueS0<TaylorModel<double> >(const std::map<std::string,TaylorModel<double> > &param) const;
template TaylorModel<Iv> ValueSet::valueS0<TaylorModel<Iv> >(const std::map<std::string,TaylorModel<Iv> > &param) const;

template <class T>
T ValueSet::get(const std::map<std::string,T> &param, const std::string &name) const
{
    typename std::map<std::string,T>::const_iterator iter = param.find(name);
    if (iter == param.end())
        throw("ValueSet::get - missing parameter!");

    return iter->second;
}

template TaylorModel<double> ValueSet::get<TaylorModel<double> >(const std::map<std::string,TaylorModel<double> > &param, const std::string &name) const;
template TaylorModel<Iv> ValueSet::get<TaylorModel<Iv> >(const std::map<std::string,TaylorModel<Iv> > &param, const std::string &name) const;

template <class T>
std::complex<T> ValueSet::s(const T &sigma, const T &omega) const
{
    if (inverse) {
        if (doBoundary()) {
            if (delta == 0.0) {
                return std::complex<T>(false,-1.0/omega);
            } else {
                return std::complex<T>(-delta,-1.0)/(omega*(delta*delta+1.0));
            }
        } else {
            if (delta == 0.0) {
                std::complex<T> ret = std::complex<T>(1.0)/std::complex<T>(sigma,omega);
                if (ret.real().bound().lower() < 0.0) {
                    std::cout << "Replacing real part of Taylor model of s^(-1) with non-negative interval!" << std::endl;
                    //std::cout << "Old:" << std::endl << ret << std::endl;
                    ret = std::complex<T>(ret.real().bound().intersect(Iv(0,INFINITY)),ret.imag());
                    //std::cout << "New:" << std::endl << ret << std::endl;
                }
                return ret;
            } else {
                return std::complex<T>(1.0)/std::complex<T>(sigma-delta*omega,omega);
            }
        }
    } else {
        if (delta == 0.0) {
            if (boundary) {
                return std::complex<T>(false,omega);
            } else {
                return std::complex<T>(sigma,omega);
            }
        } else {
            return std::complex<T>(sigma-delta*omega,omega);
        }
    }
}

template std::complex<TaylorModel<double> > ValueSet::s<TaylorModel<double> >(const TaylorModel<double> &sigma, const TaylorModel<double> &omega) const;
template std::complex<TaylorModel<Iv> > ValueSet::s<TaylorModel<Iv> >(const TaylorModel<Iv> &sigma, const TaylorModel<Iv> &omega) const;


template <class T>
std::complex<T> ValueSet::s(const T &sigma, const T &omega, size_t i, size_t n) const
{
    if (i>n)
        throw("ValueSet::s - called with i>n");

    T binom, sig, om;

    if (inverse) {
        i = n-i;
    }

    std::complex<T> ret(0.0,0.0);
    for (size_t k=0; k<=i; k++) {
        binom = T(boost::math::binomial_coefficient<double>(i,k));

        om = omega.power(k);
        sig = sigma.power(i-k);

        // actually binom*sig*om * i^k  !!!
        switch(k%4) {
        case 0:
            ret += std::complex<T>(binom*sig*om,false);
            break;
        case 1:
            ret += std::complex<T>(false,binom*sig*om);
            break;
        case 2:
            ret += std::complex<T>(-binom*sig*om,false);
            break;
        case 3:
            ret += std::complex<T>(false,-binom*sig*om);
            break;
        }
    }
    return ret;
}

template std::complex<TaylorModel<double> > ValueSet::s<TaylorModel<double> >(const TaylorModel<double> &sigma, const TaylorModel<double> &omega, size_t i, size_t n) const;
template std::complex<TaylorModel<Iv> > ValueSet::s<TaylorModel<Iv> >(const TaylorModel<Iv> &sigma, const TaylorModel<Iv> &omega, size_t i, size_t n) const;

std::ostream & operator<<(std::ostream &os, const ValueSet &vs)
{
    std::ostringstream buf;

    buf << "Value-Set:" << std::endl;
    buf << vs.ch_eq << std::endl;
    buf << "delta:     " << vs.delta << std::endl;
    buf << "boundary:  " << vs.doBoundary() << std::endl;
    buf << "complex:   " << vs.doComplex() << std::endl;
    buf << "inverse:   " << vs.doInverse() << std::endl;

    return (os << buf.str());
}
