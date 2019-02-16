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

#ifndef MONOMIAL_H
#define MONOMIAL_H

#include <iostream>
#include <string>
#include <inttypes.h>
#include <map>
#include <vector>
#include <boost/operators.hpp>

#include "iv.h"

template <class T>
class Polynomial;

template <class T>
class Monomial;

//! Stream operator for debug-output
template <class T>
std::ostream & operator<<(std::ostream &os, const Monomial<T>& m);

/*!
 * \brief The Monomial class implements a multivariate Monomial.
 *
 * A sum of Monomials makes up a multivariate Polynomial. Arithmetic operations
 * as well as functions to calculate the range of a Monomial are provided.
 *
 * The number and name of the variables can either be fixed using fixVariables
 * or can vary for each Monomial (and therefore Polynomial). The latter is much
 * more computationally complex, so it is advised to fix the variables whenever
 * possible.
 *
 * \todo Each exponent is constrained to the range uint8_t with no check to make
 * sure there is no overflow.
 */
template <class T=double>
class Monomial
        : boost::addable< Monomial<T>
        , boost::multipliable< Monomial<T>
        > >
{
public:
    //! Construct monomial with constant value coeff (exponent is zeros).
    Monomial(T coeff = 0.0);
    //! Construct monomial with one variable named var (one exponent non-zero).
    Monomial(const std::string &var);
    //! Destructor.
    ~Monomial();

    //! Return exponent
    const std::vector<uint8_t> &exp() const {return exponent;}
    //! Return coefficient
    const T &coeff() const {return coefficient;}
    //! Return variable names
    const std::vector<std::string> &vars() const {return  (variable());}
    //! Return total degree (sum of exponents)
    size_t totalDegree() const {return total_degree;}

    //! Determine value the Monomial can take if all variables are in domain.
    Iv bound(const Iv &domain) const;

    /*!
     * \brief Transform the Monomial to a new domain returning a Polynomial
     * \param scale Scale of each Monomial
     * \param shift Offset for each Monomial
     * \return Scaled and shifted Monomial (as Polynomial)
     *
     * The effect is best illustrated on a minimal example:
     * * Monomial: 1.3 a*b^2
     * * Scaled and shifted Monomial: 1.3*(shift+scale*a)*(shift+scale*b)^2
     */
    Polynomial<T> changeDomain(double scale, double shift) const;

    //! Multiply Monomial rhs with this Monomial returning my address
    Monomial &operator*= (const Monomial &rhs);
    //! Add Monomial rhs to this Monomial returning my address
    Monomial &operator+= (const Monomial &rhs);
    //! Return negation of this Monomial
    Monomial operator-() const;

    //! Comparisson of exponents to find Monomial with ideentical exponent
    bool operator<(const Monomial &rhs) const;
    //! Test whether exponents are identical
    bool operator==(const Monomial &rhs) const;

    /*!
     * \brief Fix variables resulting in much faster calculations.
     * \param var The fixed variable names that are set.
     *
     * If the number of variables is below 8 the exponent_hash is used to
     * further accelerate searching and sorting.
     */
    static void fixVariables(const std::vector<std::string> &var);
    //! Are the variables fixed?
    static bool variablesFixed() {return variable_is_fixed;}
private:
    //! Return vector of variable-names (fixed or not)
    const std::vector<std::string> &variable() const;
    //! Return map from variable-name to position in exponent (fixed or not)
    const std::map<std::string,size_t> &variableIdx() const;

    //! Calculate hash of exponent-vector
    void computeHash() const;

    //! Variables of Monomial. Empty for fixed variables.
    std::vector<std::string> var;
    //! Map from variable name to index. Empty for fixed variables.
    std::map<std::string,size_t> varidx;

    //! Monomial coefficient.
    T coefficient;
    //! Monomial exponent vector.
    std::vector<uint8_t> exponent;
    //! Hash of exponent vector
    mutable uint64_t exponent_hash;
    //! Total degree of Monomial (sum of exponent values)
    uint16_t total_degree;

    //! Merge the variable-vectors of thw Monomials
    void mergeVariables(Monomial &a, Monomial &b) const;
    //! Adapt exponent to new variable-vector and replace old variable vector
    void adjustExponent(Monomial &p, const std::vector<std::string> &var_new, const std::map<std::string,size_t> &varidx_new) const;

    //! Fixed variables. Empty for non-fixed variables.
    static std::vector<std::string> var_fixed;
    //! Map from variable name to index for fixed variables. Empty for non-fixed variables.
    static std::map<std::string,size_t> varidx_fixed;
    //! Are we using fixed variables?
    static bool variable_is_fixed;
    //! Are we using a hash to compare exponents (<=8 fixed variables)
    static bool use_exponent_hash;

    //! Stream operator for debug-output
    friend std::ostream & operator<< <> (std::ostream &os, const Monomial<T>& m);
};

#endif //MONOMIAL_H
