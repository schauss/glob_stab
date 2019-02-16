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

#ifndef BERNSTEIN_EXPONENT_H
#define BERNSTEIN_EXPONENT_H

#include "vec.hpp"
#include <map>
#include <inttypes.h>

#include <boost/thread/mutex.hpp>

/*!
 * \brief The BernsteinExponent class implements Bernstein basis polynomials of a given degree
 *
 * The Bernstein basis polynomials are used to transform a Polynomial into Bernstein form. Moreover, this class is used
 * to determine the number of Bernstein coefficients for a Bernstein polynomial of a given degree, convert an ID to an
 * exponent vector and vice versa, as well as determine the IDs of vertex coefficients.
 *
 * As the involved computations (especiall for the Bernstein basis polynomials) are quite complex the same instance
 * is used for all Bernstein instances with the same maximum degree.
 */
class BernsteinExponent
{
public:
    /*!
     * \brief Return instance of BernsteinExponent for a given maximum degree
     * \param exp maximum degree
     *
     * If no BernsteinExponent for this degree exists it is first instantiated.
     * \todo Use shared pointers!
     */
    static BernsteinExponent *getInstance(const std::vector<uint8_t> &exp);
    //! Clear all instantiated BernsteinExponents
    static void clear();

    //! Return degree of given id as vector
    Vec<uint8_t> toVec(unsigned int id) const;

    //! Return degree of given id as array
    const uint8_t *toArray(unsigned int id) const;

    //! Return ID of degree given as array
    unsigned int toId(const uint8_t exponents[]) const;

    //! Return ID of degree given as vector
    unsigned int toId(const Vec<uint8_t> &exponents) const;

    //! Maximum ID of BernsteinExponent
    unsigned int maxId() const;

    //! Check whether given ID is vertex
    bool isVertex(unsigned int id) const;

    //! Return vector of IDs of vertices
    const Vec<size_t> &vertexId() const {return vertex_id;}

    //! Return ID of linear coefficient
    const size_t *linId() const {return lin_id.data();}

    //! Number of variables
    size_t numVar() const {return max_degree.size();}

    Vec<uint8_t> max_degree; //!< Maximum degree for each dimension. \todo Make private
    Vec<size_t> step;        //!< Step size in each dimension. \todo Make private
    Vec<double> binom_vec;   //!< Vector of binomial coefficients needed to initialize Bernstein polynoial. \todo Make private

private:
    //! Construct BernsteinExponent for given maximum degree.
    BernsteinExponent(const std::vector<uint8_t> &exp);
    //! Destructor.
    ~BernsteinExponent();
    //! Instantiated BernsteinExponents stored as map with maximum degree as key
    static std::map<std::vector<uint8_t>,BernsteinExponent*> instance;
    //! Mutex to lock BernsteinExponent-map
    static boost::mutex instance_mutex;

    //! Calculate binomial coefficients
    double binomCoeffVec(const uint8_t n[], const uint8_t k[], const Vec<double> &binom, uint8_t order);

    unsigned int max_id;   //!< Maximum ID
    Vec<uint8_t> vec_id;   //!< Exponents per ID (1D for efficiency)
    Vec<size_t> vertex_id; //!< IDs of vertices
    Vec<size_t> lin_id;    //!< IDs of linear coefficients
};

#endif // BERNSTEIN_EXPONENT_H
