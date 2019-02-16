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

#ifndef MAT_HPP
#define MAT_HPP

#include "vec.hpp"

/*!
 * \brief The Mat class implements a matrix by inheriting Vec<Vec<T> >.
 *
 * The matric is stored row-wise, i.e., as Vec of rows (represented as Vec).
 * Some additional functionality is provided to return a row or column and to call functions for all members.
 *
 * \todo Get rid of this class. As inheriting from std::vector is not recomended, see, e.g.,
 * http://stackoverflow.com/questions/4353203/thou-shalt-not-inherit-from-stdvector it is probably also not a good idea
 * to inherit from a class which inherits from std::vector...
 */
template <class T>
class Mat : public Vec<Vec<T> >
{
public:
    //! Default constructor, empty matrix.
    Mat() : Vec<Vec<T> >() {}
    //! Construct matrix with given number of rows and columns setting all values to __value
    Mat(size_t rows, size_t cols, const T& __value = T()) : Vec<Vec<T> >(rows,Vec<T> (cols,__value)) {}

    //! Copy constructor. Makes use of Vec copy constructor.
    template <class U>
    Mat(const Mat<U> &m) : Vec<Vec<T> >()
    {
        for (size_t i=0; i<m.size(); i++) {
            this->push_back(Vec<T>(m[i]));
        }
    }

    /*!
     * \brief Return row with given index
     * \param index row to return (starting at 0)
     */
    Vec<T> row(size_t index)
    {
        if (index > this->size())
            throw("Out of bounds in Mat::row");

        return (*this)[index];
    }

    /*!
     * \brief Return row with given index
     * \param index column to return (starting at 0)
     *
     * Computationally more expensive than returning a row as data is stored row-wise and we need to create a new Vec
     * with one entry from each row.
     */
    Vec<T> column(size_t index)
    {
        if (index > (*this)[0].size())
            throw("Out of bounds in Mat::column");

        Vec<T> ret;
        for (size_t i=0; i < this->size(); i++)
            ret.push_back((*this)[i][index]);
        return ret;
    }
};

//! Call const function of members of elements of Mat on all elements and return Mat of results
template <class T, class U>
Mat<U> apply(const Mat<T> &mat, U (T::*mem_fn)() const)
{
    Mat<U> ret(mat.size(),mat[0].size());
    for (size_t i = 0; i<mat.size(); i++) {
        for (size_t j = 0; j<mat[0].size(); j++) {
            ret[i][j] = (mat[i][j].*mem_fn)();
        }
    }
    return ret;
}

//! Call void function of members of elements of Mat
template <class T>
void call(const Mat<T> &mat, void (T::*mem_fn)() const)
{
    for (size_t i = 0; i<mat.size(); i++) {
        for (size_t j = 0; j<mat[0].size(); j++) {
            (mat[i][j].*mem_fn)();
        }
    }
}
#endif // MAT_HPP
