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

#ifndef VEC_HPP
#define VEC_HPP

#include <cstdio>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <boost/operators.hpp>

/*!
 * \brief The Vec class implements a vector and inherits std::vector.
 *
 * In addition, element-wise arithmetic operations and comparisons, boolean checks (all, any), application of a function
 * to all elements, and several other convenient functions are implemented.
 *
 * \todo Get rid of this class. Inheriting from std::vector is not recomended, see, e.g.,
 * http://stackoverflow.com/questions/4353203/thou-shalt-not-inherit-from-stdvector
 * Switch to std::vector and algorithms (several of the ones provided here exist in the STL)
 */
template <class T>
class Vec : public std::vector<T>
        , boost::arithmetic< Vec<T> >
{
public:
    //! Default constructor, empty vector.
    Vec() : std::vector<T>() {}
    //! Construct Vec from std::vector
    Vec(const std::vector<T> &v) : std::vector<T>(v) {}
    //! Construct a vector of length 1 from one element
    Vec(const T& __value) : std::vector<T>(1,__value) {}
    //! Construct a vector of length __n with all elements set to __value
    Vec(size_t __n, const T& __value) : std::vector<T>(__n,__value) {}

    //! Copy constructor which copies all elements in v
    template <class U>
    Vec(const Vec<U> &v) : std::vector<T>()
    {
        this->reserve(v.size());
        for (size_t i=0; i<v.size(); i++)
            push_back(static_cast<T>(v[i]));
    }

    //! Construct Vec from std::vector and cast each element from type U to type T
    template <class U>
    Vec(const std::vector<U> &v) : std::vector<T>()
    {
        this->reserve(v.size());
        for (size_t i=0; i<v.size(); i++)
            push_back(static_cast<T>(v[i]));
    }

    //! Push element onto vector
    void push_back ( const T& x )
    {
        std::vector<T>::push_back(x);
    }

    //! Expand Vec by pushing each element in x
    void push_back ( const Vec<T>& x )
    {
        this->reserve(this->size()+x.size());
        for (size_t i = 0; i<x.size(); i++)
            std::vector<T>::push_back(x[i]);
    }

    //! Return Vec with opposit order of elements
    Vec<T> flip() const
    {
        Vec<T> ret;
        ret.reserve(this->size());
        for (size_t i = 1; i<=this->size(); i++)
            ret.push_back((*this)[this->size()-i]);
        return ret;
    }

    //! Comparison operations
#define compop(op)                                                  \
    Vec<bool> operator op (const Vec<T> &rhs) const                 \
    {                                                               \
        if (sizeBad(rhs))                                           \
            throw("Error in Vec:: "#op ", size missmatch!");        \
                                                                    \
        Vec<bool> ret;                                              \
        for (size_t i = 0; i<this->size(); i++)                     \
            ret.push_back((*this)[i] op rhs[rhs.size()!=1?i:0]);    \
                                                                    \
        return ret;                                                 \
    }                                                               \
                                                                    \
    Vec<bool> operator op (const T &rhs) const                      \
    {                                                               \
        Vec<T> rhsv(1,rhs);                                         \
        return (*this op rhsv);                                     \
    }

    compop(<)
    compop(<=)
    compop(>)
    compop(>=)
    compop(==)
    compop(!=)
#undef compop

    //! Arithmetic operations
#define operatorArith(op)                                           \
    Vec<T> &operator op ## = (const Vec<T> &rhs)                    \
    {                                                               \
        if (sizeBad(rhs)) {                                         \
            throw("Error in Vec:: "#op "=, size missmatch!");       \
        }                                                           \
                                                                    \
        Vec<T> ret;                                                 \
        for (size_t i = 0; i<this->size(); i++)                     \
            (*this)[i] op ## = rhs[rhs.size()!=1?i:0];              \
                                                                    \
        return *this;                                               \
    }

    operatorArith(+)
    operatorArith(-)
    operatorArith(*)
    operatorArith(/)
#undef operatorArith

    //! Negation of all elements
    Vec<T> operator- () const
    {
        Vec<T> ret;
        for (size_t i = 0; i<this->size(); i++)
            ret.push_back(-(*this)[i]);

        return ret;
    }

    //! Logical invers of all elements
    Vec<bool> operator! () const
    {
        Vec<T> ret;
        for (size_t i = 0; i<this->size(); i++)
            ret.push_back(!(*this)[i]);

        return ret;
    }

    //! Assignment operator
    Vec<T> &operator=(const Vec<T> &rhs)
    {
        this->clear();
        for (size_t i = 0; i<rhs.size(); i++)
             push_back(rhs[i]);
        return *this;
    }

    //! Assign single value => Vec of length one
    Vec<T> &operator=(const T &rhs)
    {
        this->clear();
        push_back(rhs);
        return *this;
    }

    //! Returns true if all elements are true
    bool all(int *first=0) const
    {
        bool ret = true;
        for (size_t i = 0; i<this->size(); i++) {
            if ((*this)[i])
                continue;
            else {
                ret = false;
                if (first)
                    *first = i;
                break;
            }
        }
        return ret;
    }

    //! Returns true if any element is true
    bool any(int *first=0) const
    {
        bool ret = false;
        for (size_t i = 0; i<this->size(); i++) {
            if ((*this)[i]) {
                ret = true;
                if (first)
                    *first = i;
                break;
            }
        }
        return ret;
    }

    //! Return iterator to minimum element (and index of minimum element if a non-zero pointer is provided)
    typename Vec<T>::iterator min(int *index=0)
    {
        // typename Vec<T>::iterator ret = std::min_element(this->begin(),this->end());
        // for double this makes shure nan's are ignored...
        typename Vec<T>::iterator ret = this->begin();
        for (typename Vec<T>::iterator i = this->begin()+1; i<this->end(); i++) {
            if ((*i < *ret) || std::isnan(*ret) )
                ret = i;
        }
        if (index)
            *index = ret - this->begin();
        return ret;
    }

    //! Return value of minimum element (and index of minimum element if a non-zero pointer is provided)
    T minVal(int *index=0)
    {
        typename Vec<T>::iterator ret = min(index);
        return *ret;
    }

    //! Return iterator to maximum element (and index of minimum element if a non-zero pointer is provided)
    typename Vec<T>::iterator max(int *index=0)
    {
        // typename Vec<T>::iterator ret = std::max_element(this->begin(),this->end());
        // for double this makes shure nan's are ignored...
        typename Vec<T>::iterator ret = this->begin();
        for (typename Vec<T>::iterator i = this->begin()+1; i<this->end(); i++) {
            if ((*i > *ret) || std::isnan(*ret) )
                ret = i;
        }
        if (index)
            *index = ret - this->begin();
        return ret;
    }

    //! Return value of maximum element (and index of minimum element if a non-zero pointer is provided)
    T maxVal(int *index=0)
    {
        typename Vec<T>::iterator ret = max(index);
        return *ret;
    }

    //! Sum of all elements of Vec
    T sum() const
    {
        T ret = 0.0;
        for (size_t i = 0; i<this->size(); i++)
            ret += (*this)[i];

        return ret;
    }

    //! Cumulative sum of elements of Vec (same length as Vec)
    Vec<T> cumSum() const
    {
        Vec<T> ret = *this;
        for (size_t i = 1; i<this->size(); i++)
            ret[i] += ret[i-1];

        return ret;
    }

    //! Product of all elements of Vec
    T prod() const
    {
        T ret = 1.0;
        for (size_t i = 0; i<this->size(); i++)
            ret *= (*this)[i];

        return ret;
    }

    //! Cumulative product of elements of Vec (same length as Vec)
    Vec<T> cumProd() const
    {
        Vec<T> ret = *this;
        for (size_t i = 1; i<this->size(); i++)
            ret[i] *= ret[i-1];

        return ret;
    }

    //! Check whether size of Vec and rhs are equal
    template <class U>
    bool sizeBad(const Vec<U> &rhs) const
    {
        return ( (rhs.size() != this->size()) && (rhs.size() != 1) );
    }
};

//! Call function of members of elements of Vec on all elements and return Vec of results
template <class T, class U>
Vec<U> apply(Vec<T> vec, U (T::*mem_fn)())
{
    Vec<U> ret;
    for (size_t i = 0; i<vec.size(); i++) {
        ret.push_back( (vec[i].*mem_fn)() );
    }
    return ret;
}

/*
template <class T, class U>
Vec<U> apply(Vec<T> &vec, U (T::*mem_fn)())
{
    Vec<U> ret;
    for (size_t i = 0; i<vec.size(); i++) {
        ret.push_back( (vec[i].*mem_fn)() );
    }
    return ret;
}
*/

//! Call const function of members of elements of Vec on all elements and return Vec of results
template <class T, class U>
Vec<U> apply(const Vec<T> &vec, U (T::*mem_fn)() const)
{
    Vec<U> ret;
    for (size_t i = 0; i<vec.size(); i++) {
        ret.push_back( (vec[i].*mem_fn)() );
    }
    return ret;
}

//! Call function for all elements of Vec and return Vec of results
template <class T, class U>
Vec<U> apply(const Vec<T> &vec, U (*fn)(T))
{
    Vec<U> ret;
    for (size_t i = 0; i<vec.size(); i++) {
        ret.push_back( ( *fn)(vec[i]) );
    }
    return ret;
}

//! Call void function of members of elements of Vec
template <class T>
void call(const Vec<T> &vec, void (T::*mem_fn)() const)
{
    for (size_t i = 0; i<vec.size(); i++) {
        (vec[i].*mem_fn)();
    }
}

//! Return Vec of given member-variable of all elements of Vec
template <class T, class U>
Vec<U> getVec(const Vec<T> &vec, U (T::*mem_val) )
{
    Vec<U> ret;
    for (size_t i = 0; i<vec.size(); i++)
        ret.push_back(vec[i].*mem_val);

    return ret;
}

//! Set member-variable of all elements of Vec to values in val
template <class T, class U, class V>
void setVec(Vec<T> &vec, U (T::*mem_val), Vec<V> val)
{
    if (vec.sizeBad(val))
         throw("Error in Vec::set, size missmatch!");

    for (size_t i = 0; i<vec.size(); i++)
        vec[i].*mem_val = val[i];
}

//! Return Vec of sine-values of elements of vec
template <class T>
Vec<T> sin(const Vec<T> &vec)
{
    Vec<T> ret;
    for (size_t i = 0; i<vec.size(); i++) {
        ret.push_back(sin(vec[i]));
    }
    return ret;
}

//! Return Vec of cosine-values of elements of vec
template <class T>
Vec<T> cos(const Vec<T> &vec)
{
    Vec<T> ret;
    for (size_t i = 0; i<vec.size(); i++) {
        ret.push_back(cos(vec[i]));
    }
    return ret;
}

//! Return Vec of exponential function of elements of vec
template <class T>
Vec<T> exp(const Vec<T> &vec)
{
    Vec<T> ret;
    for (size_t i = 0; i<vec.size(); i++) {
        ret.push_back(exp(vec[i]));
    }
    return ret;
}

//! Stram-operator for output
template<class T>
std::ostream & operator<<(std::ostream &os, const Vec<T>& vec)
{
    if (vec.empty())
        return os;

    os << "[ ";

    for (size_t i = 0; i<vec.size()-1; i++)
        os << vec[i] << "; ";

    os << vec[vec.size()-1] << " ]";

    return os;
}

#endif // VEC_HPP
