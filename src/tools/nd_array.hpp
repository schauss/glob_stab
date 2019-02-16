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

#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

template <class T> class ndArrayView;

/*!
 * \brief The ndArray template class implements a n-dimensional array.
 *
 * The number of dimensions n and size in each dimension do not have to be known at compile time and are specified on
 * initialization.
 *
 * A versatile [] operator is provided which allows to access single elements, n-1-dimensional arrays, or ranges of
 * elements.
 */
template <class T>
class ndArray
{
public:
    //! Create empty ndArray
    ndArray() : dim(0)
    {
    }

    //! Create ndArray with size of dimensions given in vector dimension_size
    ndArray(std::vector<size_t> &dimension_size) :
        dim(dimension_size.size()), dim_fixed(0), dim_size(dimension_size),
        dim_start(std::vector<size_t>(dim,0)), dim_stride(std::vector<size_t>(dim,1)), dim_end(dim_size)
    {
        size_t data_size = 1;
        for (size_t i=0; i<dim; i++) {
            if (dim_size[i] == 0)
                throw("ndArray::ndArray() - dimension has size 0!");
            data_size *= dim_size[i];
        }
        data = boost::shared_ptr<std::vector<T> >(new std::vector<T>(data_size));
    }

    //! Set all values of the ndArray to val. See ndArrayView::operator=
    void operator=(const T& val)
    {
        ndArrayView<T>(*this) = val;
    }

    //! Creates an ndArrayView from an ndArray.
    operator T() const
    {
        return ndArrayView<T>(*this);
    }

    //! Reduces the dimension by one, selecting idx in the first "free" dimension. See ndArrayView::operator[]
    ndArrayView<T> operator[](size_t idx) const
    {
        return ndArrayView<T>(*this)[idx];
    }

    /*!
     * \brief Keeps the number of dimensions constant, returning a n-dimensional array where each dimension may be
     * smaller than of the original array. See ndArrayView::operator()
     * \param start_idx Start index for each dimension
     * \param end_idx End index for each dimension
     */
    ndArrayView<T> operator()(const std::vector<size_t> &start_idx, const std::vector<size_t> &end_idx) const
    {
        return ndArrayView<T>(*this)(start_idx,end_idx);
    }

    //! Returns true if all values in the ndArray have the value val. See ndArrayView::operator==
    bool operator==(const T& val) const
    {
        return (ndArrayView<T>(*this) == val);
    }

    //! Returns true if any value in the ndArray has the value val. See ndArrayView::any()
    bool any(const T& val) const
    {
        return (ndArrayView<T>(*this).any(val));
    }

    //! Returns all indices for which the ndArray has the value val (as vector of vectors). See ndArrayView::find()
    std::vector<std::vector<size_t> > find(const T& val) const
    {
        return (ndArrayView<T>(*this).find(val));
    }

    //! See ndArrayView::boundingBox()
    bool boundingBox(const T& val, std::vector<size_t> &start, std::vector<size_t> &end) const
    {
        return (ndArrayView<T>(*this).boundingBox(val,start,end));
    }

    //! Output information an ndArray to std::cout
    void print() const
    {
        std::cout << "nd-Array:" << std::endl;
        std::cout << "dim = " << dim << ", dim_fixed = " << dim_fixed << std::endl;
        for(size_t i=0; i<dim; i++) {
            std::cout << "dim " << i <<": size = " << dim_size[i] << ", slice = [" <<
                         dim_start[i] << ", " << dim_stride[i] << ", " << dim_end[i] <<
                         "]" << std::endl;
        }
    }

    friend class ndArrayView<T>;
private:
    size_t dim;                              //!< Number of dimensions
    size_t dim_fixed;                        //!< Number of fixed dimensions
    std::vector<size_t> dim_size;            //!< Size of each dimension

    std::vector<size_t> dim_start;           //!< Starting index for each dimension
    std::vector<size_t> dim_stride;          //!< Stride for each dimension
    std::vector<size_t> dim_end;             //!< End-idx for each dimension (one past last element)

    boost::shared_ptr<std::vector<T> > data; //!< Actual data as shared ptr, i.e., views of one array can use same data
};

/*!
 * \brief The ndArrayView template class implements a view onto an ndArray.
 *
 * Much of the functionality of ndArray is actually implemented in this class.
 *
 * A versatile [] operator is provided which allows to access single elements, n-1-dimensional arrays, or ranges of
 * elements.
 */
template <class T>
class ndArrayView
{
public:
    //! Construct an ndArrayView onto an ndArray
    ndArrayView(const ndArray<T> &parent) : array(parent)
    {
        if (array.dim==0)
            throw("ndArrayView::ndArrayView(parent) - parent empty!\n");
    }

    /*!
     * \brief Reduces the dimension by one by selecting the index in the first "free" dimension
     * \param idx Index value in the first "free" dimension
     */
    ndArrayView<T>& operator[](size_t idx)
    {
        if (array.dim_fixed==array.dim)
            throw("ndArrayStride::ndArrayStride(parent,idx) - all dimensions fixed!");

        array.dim_start[array.dim_fixed] += idx;
        array.dim_end[array.dim_fixed] = array.dim_start[array.dim_fixed]+1;
        array.dim_fixed++;

        return *this;
    }

    /*!
     * \brief Keeps the number of dimensions constant but manipulates the values dim_start and dim_end for all "free"
     * dimensions.
     * \param start_idx Start index for each dimension (relative to current view)
     * \param end_idx End index for each dimension (relative to current view)
     *
     * The resulting array has the same number of "free" dimensions but may have a reduced length in each dimension.
     */
    ndArrayView<T>& operator()(const std::vector<size_t> &start_idx, const std::vector<size_t> &end_idx)
    {
        if ( (start_idx.size() != array.dim-array.dim_fixed) || (end_idx.size() != array.dim-array.dim_fixed) )
            throw("ndArrayView<T> operator(std::vector,std::vector) called with wrong size vector!");

        for (size_t i=array.dim_fixed; i<array.dim; i++) {
            if (start_idx[i] >= end_idx[i])
                throw("ndArrayView<T> operator(std::vector,std::vector) - start_idx>=end_idx!");

            if (array.dim_start[i]+end_idx[i]*array.dim_stride[i] > array.dim_end[i])
                throw("ndArrayView<T> operator(std::vector,std::vector) - out of bounds!");

            array.dim_end[i] = array.dim_start[i]+end_idx[i]*array.dim_stride[i];
            array.dim_start[i] += start_idx[i]*array.dim_stride[i];
        }

        return *this;
    }

    //! Sets all values in current view to val
    void operator=(const T& val)
    {
        if (array.dim_fixed<array.dim) {
            for (size_t i = 0; i < array.dim_end[array.dim_fixed]-array.dim_start[array.dim_fixed]; i+=array.dim_stride[array.dim_fixed]) {
                array[i] = val;
            }
        }
        else
            (*(array.data.get()))[computeIdx()] = val;
    }

    //! Returns true if all values in current view compare equal to val.
    bool operator==(const T& val) const
    {
        if (array.dim_fixed<array.dim) {
            for (size_t i = 0; i < array.dim_end[array.dim_fixed]-array.dim_start[array.dim_fixed]; i+=array.dim_stride[array.dim_fixed]) {
                if (!(array[i] == val))
                    return false;
            }
            return true;
        }
        else
            return (*(array.data.get()))[computeIdx()] == val;
    }

    //! Returns true if any value in current view compares equal to val
    bool any(const T& val) const
    {
        if (array.dim_fixed<array.dim) {
            for (size_t i = 0; i < array.dim_end[array.dim_fixed]-array.dim_start[array.dim_fixed]; i+=array.dim_stride[array.dim_fixed]) {
                if (array[i].any(val))
                    return true;
            }
            return false;
        }
        else
            return (*(array.data.get()))[computeIdx()] == val;
    }


    //! Returns vector of vectors of indices of elements which are equal to val
    std::vector<std::vector<size_t> > find(const T& val) const
    {
        if (array.dim_fixed<array.dim) {
            std::vector<std::vector<size_t> > ret;
            for (size_t i = 0; i < array.dim_end[array.dim_fixed]-array.dim_start[array.dim_fixed]; i+=array.dim_stride[array.dim_fixed]) {
                std::vector<std::vector<size_t> > tmp = array[i].find(val);
                for (size_t j=0; j<tmp.size(); j++)
                    ret.push_back(tmp[j]);
            }
            return ret;
        }
        else {
            if ( (*(array.data.get()))[computeIdx()] == val )
                return std::vector<std::vector<size_t> >(1,array.dim_start);
            else
                return std::vector<std::vector<size_t> >();
        }
    }

    //! Return true if any elements equal to val, false if not, set bounding box of indexes for which elements are equal to val
    bool boundingBox(const T& val, std::vector<size_t> &start, std::vector<size_t> &end) const
    {
        std::vector<std::vector<size_t> >  eq_vec = find(val);
        if (eq_vec.empty()) {
            return false;
        }
        else {
            start = array.dim_end;
            end = array.dim_start;
            for (size_t i=0; i<eq_vec.size(); i++)
            {
                for (size_t dim=0; dim<array.dim; dim++) {
                    if (start[dim] > eq_vec[i][dim])
                        start[dim] = eq_vec[i][dim];
                    if (end[dim]   < eq_vec[i][dim]+1)
                        end[dim]   = eq_vec[i][dim]+1;
                }
            }

            return true;
        }
    }

    //! Return one element of array (all dimensions must be fixed).
    operator T() const
    {
        if (array.dim_fixed<array.dim)
            throw("ndArray::operator T() - conversion to type T not possible, not all dimensions fixed!");

        return (*(array.data.get()))[computeIdx()];
    }

private:
    //! Calculate idx of value in data which corresponds to fixed array view
    size_t computeIdx() const
    {
        size_t idx=0;
        size_t step=1;

        for(size_t i=0; i<array.dim; i++) {
            idx += array.dim_start[i]*step;
            step *= array.dim_size[i];
        }
        return idx;
    }

    ndArray<T> array;
};

#endif // NDARRAY_HPP


