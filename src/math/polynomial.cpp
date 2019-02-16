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

#include "polynomial.h"
#include "bernstein.h"
#include <algorithm>

template <class T>
bool Polynomial<T>::do_bernstein_bound = true;

template <class T>
Polynomial<T>::Polynomial()
{
}

template <class T>
Polynomial<T>::Polynomial(const std::string& var)
{
    Monomial<T> m(var);
    monomial.insert(m);
}

template <class T>
Polynomial<T>::Polynomial(const char var[])
{
    std::string str = var;
    Monomial<T> m(str);
    monomial.insert(m);
}

template <class T>
Polynomial<T>::Polynomial(T val)
{
    if (val != 0.0) {
        Monomial<T> m(val);
        monomial.insert(m);
    }
}

template <class T>
Polynomial<T>::~Polynomial()
{

}

template <class T>
T Polynomial<T>::constant() const
{
    if ( monomial.empty() )
        return 0.0;
    else {
        typename std::set<Monomial<T> >::const_iterator it = monomial.find(Monomial<T>(0.0));
        if (it == monomial.end())
            return 0.0;
        else
            return it->coeff();
    }
}

template <class T>
const Iv &Polynomial<T>::bound(const Iv &domain, bool do_bernstein) const
{
    if (poly_bound.isNaN() || !(domain == bound_domain) || (did_bernstein != do_bernstein)) {
        bound_domain = domain;
        did_bernstein = do_bernstein;
        poly_bound = Iv(0,0);
        for (typename std::set<Monomial<T> >::const_iterator it = monomial.begin(); it!=monomial.end(); it++) {
            poly_bound += it->bound(domain);
        }
        if (do_bernstein) {
            Iv b_bound;
            if (!(domain == Iv(0,1))) {
                Bernstein<T> b(changeDomain(domain,Iv(0,1)));
                b_bound = b.bound();
            } else {
                Bernstein<T> b(*this);
                b_bound = b.bound();
            }
            //std::cout << "Polynomial::bound():" << std::endl;
            //std::cout << "iv-bound = " << poly_bound << std::endl;
            //std::cout << "b-bound  = " << b_bound << std::endl;
            poly_bound = poly_bound.intersect(b_bound);
            //std::cout << "bound    = " << b_bound << std::endl << std::endl;
        }
    }
    return poly_bound;
}

template <class T>
Iv Polynomial<T>::cutoff(size_t total_degree, const Iv &domain)
{
    Iv ret(false);

    if (monomial.empty())
        return ret;

    for (typename std::set<Monomial<T> >::iterator it = monomial.begin(); it!=monomial.end();) {
        if (it->totalDegree() > total_degree) {
            ret += it->bound(domain);
            monomial.erase(it++);
        }
        else
            ++it;
    }

    return ret;
}

template <class T>
Polynomial<T> Polynomial<T>::changeDomain(const Iv &old_dom, const Iv &new_dom) const
{
    double scale = old_dom.width()/new_dom.width();
    double shift = old_dom.lower()- scale*new_dom.lower();

    Polynomial<T> ret;

    for (typename std::set<Monomial<T> >::const_iterator it = monomial.begin(); it!=monomial.end(); it++)
        ret += it->changeDomain(scale,shift);

    return ret;
}

template <class T>
std::vector<uint8_t> Polynomial<T>::maxOrders() const
{
    if (!Monomial<T>::variablesFixed())
        throw("Polynomial::maxOrders - Maximum orders can only be determined for polynomial with fixed variables!");

    if (monomial.empty())
        return std::vector<uint8_t>(Monomial<T>().exp().size(),0);

    std::vector<uint8_t> ret(monomial.begin()->exp().size(),0);
    for (typename std::set<Monomial<T> >::const_iterator it = monomial.begin(); it!=monomial.end(); it++) {
        for (size_t e=0; e<it->exp().size(); e++) {
            if (it->exp()[e] > ret[e])
                ret[e] = it->exp()[e];
        }
    }

    return ret;
}

/*!
 * Searches for Monomials in rhs with same exponent as Monomials in this
 * Polynomial. If this search is successfull the coefficients are added. If not,
 * the new Monomial from rhs is added to this Polynomial.
 *
 * \todo Why do we first delete and then reinsert the monomial if it is found?
 */
template <class T>
Polynomial<T> &Polynomial<T>::operator+= (const Polynomial<T> &rhs)
{
    poly_bound = Iv();
    for (typename std::set<Monomial<T> >::const_iterator it_rhs = rhs.monomial.begin(); it_rhs!=rhs.monomial.end(); it_rhs++) {
        std::pair<typename std::set<Monomial<T> >::iterator,bool> add = monomial.insert(*it_rhs);
        if (!(add.second)) {
            Monomial<T> m = *(add.first) + *it_rhs;
            monomial.erase(add.first);
            if (m.coeff() != 0.0)
                monomial.insert(m);
        }
    }

    return *this;
}

/*!
 * Adds the negated rhs Polynomial
 */
template <class T>
Polynomial<T> &Polynomial<T>::operator-= (const Polynomial<T> &rhs)
{
    return (*this += -rhs);
}

/*!
 * Store my set of Monomials lhs and clear my set of Monomials. Then multiply
 * each Monomial in lhs with each Monomial in rhs. Searches for resulting
 * Monomial in my set of Monomials. If this search is successfull the
 * coefficients are added. If not, the new Monomial is added to this Polynomial.
 *
 * \todo Why do we first delete and then reinsert the monomial if it is found?
 */
template <class T>
Polynomial<T> &Polynomial<T>::operator*= (const Polynomial<T> &rhs)
{
    poly_bound = Iv();
    std::set<Monomial<T> > lhs_monomial(monomial);
    monomial = std::set<Monomial<T> >();

    for (typename std::set<Monomial<T> >::iterator it = lhs_monomial.begin(); it!=lhs_monomial.end(); it++) {
        for (typename std::set<Monomial<T> >::const_iterator it_rhs = rhs.monomial.begin(); it_rhs!=rhs.monomial.end(); it_rhs++) {
            Monomial<T> m = (*it) * (*it_rhs);
            std::pair<typename std::set<Monomial<T> >::iterator,bool> add = monomial.insert(m);
            if (!(add.second)) {
                m += *(add.first);
                monomial.erase(add.first);
                if (m.coeff() != 0.0)
                    monomial.insert(m);
            }
        }
    }

    return *this;
}

template <class T>
Polynomial<T> &Polynomial<T>::operator+= (double rhs)
{
    *this += Polynomial<T>(rhs);
    return *this;
}

template <class T>
Polynomial<T> &Polynomial<T>::operator-= (double rhs)
{
    *this -= Polynomial<T>(rhs);
    return *this;
}

template <class T>
Polynomial<T> &Polynomial<T>::operator*= (double rhs)
{
    *this *= Polynomial<T>(rhs);
    return *this;
}

template <class T>
Polynomial<T> Polynomial<T>::operator-() const
{
    Polynomial<T> ret;
    for (typename std::set<Monomial<T> >::const_iterator it = monomial.begin(); it!=monomial.end(); it++) {
        ret.monomial.insert(-(*it));
    }
    ret.poly_bound = -poly_bound;
    ret.bound_domain = bound_domain;

    return ret;
}

template <class T>
Polynomial<T> Polynomial<T>::power(unsigned int i) const
{
    if (i==0)
        return Polynomial(1.0);
    else if (i==1)
        return *this;
    else {
        if (isEmpty())
            return Polynomial(0.0);

        Polynomial ret = *this;
        for (size_t k=1; k<i;k++)
            ret *= *this;
        return ret;
    }
}

template class Polynomial<double>;
template class Polynomial<Iv>;

template <class T>
std::ostream & operator<<(std::ostream &os, const Polynomial<T>& p)
{
    if (p.monomial.empty())
        os << Monomial<T>(0.0);
    else {
        for (typename std::set<Monomial<T> >::const_iterator it = p.monomial.begin(); it!=p.monomial.end(); it++)
            os << *it;
    }

    return os;
}
template std::ostream & operator<<(std::ostream &os, const Polynomial<double>& p);
template std::ostream & operator<<(std::ostream &os, const Polynomial<Iv>& p);
