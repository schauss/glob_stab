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

#include "monomial.h"
#include "polynomial.h"

#include <algorithm>
#include <cmath>

template <class T>
std::vector<std::string> Monomial<T>::var_fixed;

template <class T>
std::map<std::string,size_t> Monomial<T>::varidx_fixed;

template <class T>
bool Monomial<T>::variable_is_fixed = false;

template <class T>
bool Monomial<T>::use_exponent_hash = false;

template <class T>
Monomial<T>::Monomial(T coeff)
    : coefficient(coeff), exponent_hash(0), total_degree(0)
{
    if (variable_is_fixed)
        exponent = std::vector<uint8_t>(variable().size(), 0);
}

template <class T>
Monomial<T>::Monomial(const std::string &var)
    : coefficient(1.0), total_degree(1)
{
    if (variable_is_fixed) {
        std::map<std::string,size_t>::const_iterator it = variableIdx().find(var);
        if (it == variableIdx().end())
            throw("Monomial::Monomial - could not find variable in varidx_fixed!");
        exponent = std::vector<uint8_t>(variable().size(), 0);
        exponent[it->second] = 1;
        if (use_exponent_hash)
            computeHash();
    } else {
        exponent.push_back(1);
        this->var.push_back(var);
        this->varidx[var] = 0;
    }
}

template <class T>
void Monomial<T>::fixVariables(const std::vector<std::string> &var)
{
    variable_is_fixed = true;
    use_exponent_hash = (var.size()<=8);
    var_fixed = var;
    for (size_t i=0; i< var.size(); i++)
        varidx_fixed[var[i]] = i;
}

template <class T>
Iv Monomial<T>::bound(const Iv &domain) const
{
    if (total_degree == 0)
        return Iv(coefficient);

    if (domain == Iv(0,1)) {
        return coefficient * Iv(0.0,1.0);
    } else if (domain == Iv(-1,1)) {

        // check for odd exponents
        bool has_odd = false;
        for (size_t i=0; i<exponent.size();i++) {
            if (exponent[i]%2 == 1) {
                has_odd = true;
                break;
            }
        }

        if (has_odd)
            return coefficient * Iv(-1.0,1.0);
        else
            return coefficient * Iv(0.0,1.0);

    } else {
        throw("Monomial::bound - not defined for domain other than [0,1] or [-1,1]");
    }
}

// for doubles we can use a special "quick solution" for domain = Iv(0,1)
template <>
Iv Monomial<double>::bound(const Iv &domain) const
{
    if (total_degree == 0)
        return Iv(coefficient);

    if (domain == Iv(0,1)) {
        if (coefficient >= 0.0)
            return Iv(0.0,coefficient);
        else
            return Iv(coefficient,0.0);
    } else if (domain == Iv(-1,1)) {

        // check for odd exponents
        bool has_odd = false;
        for (size_t i=0; i<exponent.size();i++) {
            if (exponent[i]%2 == 1) {
                has_odd = true;
                break;
            }
        }

        if (has_odd)
            return coefficient * Iv(-1.0,1.0);
        else
            return coefficient * Iv(0.0,1.0);

    } else {
        throw("Monomial::bound - not defined for domain other than [0,1] or [-1,1]");
    }
}

template <class T>
Polynomial<T> Monomial<T>::changeDomain(double scale, double shift) const
{
    Polynomial<T> ret(coefficient);
    for (size_t i=0; i<exponent.size(); i++) {
        if (exponent[i] == 0)
            continue;
        ret *= (Polynomial<T>(shift) + Polynomial<T>(scale)*Polynomial<T>(variable()[i])).power(exponent[i]);
    }
    return ret;
}

/*!
 * The coefficients are multiplied and the exponents are added.
 */
template <class T>
Monomial<T> &Monomial<T>::operator*= (const Monomial<T> &rhs)
{
    if (variable_is_fixed) {
        coefficient *= rhs.coefficient;
        for (size_t i=0; i<exponent.size();i++)
            exponent[i] += rhs.exponent[i];
        total_degree += rhs.total_degree;
        if (use_exponent_hash)
            computeHash();
    } else {
        Monomial<T> r(rhs);
        mergeVariables(*this,r);

        coefficient *= r.coefficient;
        for (size_t i=0; i<exponent.size();i++)
            exponent[i] += r.exponent[i];
        total_degree += rhs.total_degree;
    }
    return *this;
}

/*!
 * The coefficients are added it the exponents are identical, otherwise an
 * Exception is raised.
 */
template <class T>
Monomial<T> &Monomial<T>::operator+= (const Monomial<T> &rhs)
{
    if (!(*this == rhs))
        throw("Monomial::+ - exponent missmatch!");

    coefficient += rhs.coefficient;
    return *this;
}

template <class T>
Monomial<T> Monomial<T>::operator-() const
{
    Monomial<T> ret(*this);
    ret.coefficient = -coefficient;

    return ret;
}

template <class T>
bool Monomial<T>::operator<(const Monomial<T> &rhs) const
{
    if (use_exponent_hash)
        return (exponent_hash < rhs.exponent_hash);

    if (total_degree < rhs.total_degree)
        return true;
    else if (total_degree > rhs.total_degree)
        return false;

    if (!variable_is_fixed) {
        if (exponent.size() < rhs.exponent.size())
            return true;
        else if (exponent.size() > rhs.exponent.size())
            return false;

        for (size_t i=0; i<variable().size(); i++) {
            if (variable()[i] < rhs.variable()[i])
                return true;
            else if (variable()[i] > rhs.variable()[i])
                return false;
        }
    }

    for (int i=exponent.size()-1; i>=0; i--) {
        if (exponent[(size_t)i] < rhs.exponent[(size_t)i])
            return true;
        else if (exponent[(size_t)i] > rhs.exponent[(size_t)i])
            return false;
    }

    // both the same!
    return false;
}

template <class T>
bool Monomial<T>::operator==(const Monomial<T> &rhs) const
{
    if (!variable_is_fixed && (exponent.size() != rhs.exponent.size()) )
            return false;

    for (size_t i=0; i<exponent.size(); i++) {
        if (!variable_is_fixed && (variable()[i] != rhs.variable()[i]) )
            return false;

        if (exponent[i] != rhs.exponent[i])
            return false;
    }

    // both the same!
    return true;
}

template <class T>
Monomial<T>::~Monomial()
{

}

template <class T>
const std::vector<std::string> &Monomial<T>::variable() const
{
    if (variable_is_fixed)
        return var_fixed;
    else
        return var;
}

template <class T>
const std::map<std::string,size_t> &Monomial<T>::variableIdx() const
{
    if (variable_is_fixed)
        return varidx_fixed;
    else
        return varidx;
}

/*!
 * Calculates one uint64_t from up to eight uint8_t. Assumes that the exponent-
 * size is at most eight (this is not checked).
 *
 * \todo Could we simply use a union? If this works it would be more efficient!
 */
template <class T>
void  Monomial<T>::computeHash() const
{
    exponent_hash = exponent[0];
    for (size_t i=1; i<exponent.size(); i++) {
        exponent_hash = (exponent_hash<<8)+exponent[i];
    }
}

/*!
 * Computationally quite expensive. Not necessary for fixed variables!
 */
template <class T>
void  Monomial<T>::mergeVariables(Monomial<T> &a,  Monomial<T> &b) const
{
    // fast exit if both variable-vectors are identical
    if (a.variable().size() == b.variable().size()) {
        bool done = true;
        for (size_t i=0; i<a.variable().size(); i++) {
            if (a.variable()[i] != b.variable()[i]) {
                done = false;
                break;
            }
        }
        if (done)
            return;
    }

    // otherwise merge two variable-vectors
    std::vector<std::string> variable_merged(a.variable().size()+b.variable().size(),"");
    std::merge(a.variable().begin(),a.variable().end(),b.variable().begin(),b.variable().end(),variable_merged.begin());
    std::vector<std::string>::iterator it = std::unique(variable_merged.begin(),variable_merged.end());
    variable_merged.resize(it - variable_merged.begin());

    std::map<std::string,size_t> varidx_merged;
    for (size_t i=0; i<variable_merged.size(); i++) {
        varidx_merged[variable_merged[i]] = i;
    }

    // adjust the exponent-vectors accordingly and replace variable-vectors
    adjustExponent(a, variable_merged, varidx_merged);
    adjustExponent(b, variable_merged, varidx_merged);
}

/*!
 * Computationally quite expensive. Not necessary for fixed variables!
 */
template <class T>
void  Monomial<T>::adjustExponent(Monomial<T> &p, const std::vector<std::string> &var_new, const std::map<std::string,size_t> &varidx_new) const
{
    std::vector<uint8_t> exp(var_new.size(),0);
    for (size_t i = 0; i < p.variable().size(); i++) {
        std::map<std::string,size_t>::const_iterator it = varidx_new.find(p.variable()[i]);
        if (it == varidx_new.end())
            throw("Monomial::adjustExponent - could not find variable in varidx_new!");
        exp[it->second] = p.exponent[i];
    }
    p.exponent = exp;

    p.var = var_new;
    p.varidx = varidx_new;
}

template class Monomial<double>;
template class Monomial<Iv>;

template <class T>
std::ostream & operator<<(std::ostream &os, const Monomial<T>& m)
{

    os << (m.coefficient>=0.0?" + ":" - ") << fabs(m.coefficient);

    for (size_t i = 0; i<m.exponent.size(); i++) {
        if (m.exponent[i] > 0) {
            os << " " << (m.variable()[i]);
            if (m.exponent[i] > 1) {
                os << "^" << (int)m.exponent[i];
            }
        }
    }

    return os;
}
template std::ostream & operator<<(std::ostream &os, const Monomial<double>& m);

template <>
std::ostream & operator<<(std::ostream &os, const Monomial<Iv>& m)
{

    os << " + " << m.coefficient;

    for (size_t i = 0; i<m.exponent.size(); i++) {
        if (m.exponent[i] > 0) {
            os << " " << (m.variable()[i]);
            if (m.exponent[i] > 1) {
                os << "^" << (int)m.exponent[i];
            }
        }
    }

    return os;
}
template std::ostream & operator<<(std::ostream &os, const Monomial<Iv>& m);
