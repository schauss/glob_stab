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

#ifndef EQUATION_H
#define EQUATION_H

#include "vec.hpp"
#include "mat.hpp"
#include <complex>
#include <map>
#include <string>
#include <ostream>

/*!
 * \brief The class Equation evaluates a mathematical expression.
 *
 * The class makes use of ParseEquation to evaluate the equation on calls to value() or valueMat().
 */
class Equation
{
public:
    /*!
     * \brief Constructor which does not perform any evaluation. Simply sets the equation string and name of standard
     * return variable.
     * \param equation Set member eq to equation
     * \param variable Set member var to variable
     */
    Equation(const std::string &equation, const std::string &variable);
    //! Destructor.
    virtual ~Equation() {}

    /*!
     * \brief Evaluate equation and return result as Vec
     * \param parameters Parameter values to use for evaluation
     * \param variable Variable to return, if not specified then the variable specified in the constructor and stored in
     * var is returned
     *
     * If the same variables exist in parameters and in fixed_vals (set using set()), then the values from fixed_vals
     * are used.
     */
    template <class T>
    Vec<T> value(const std::map<std::string,T> &parameters, const std::string &variable = std::string()) const;

    /*!
     * \brief Evaluate equation and return result as Mat
     * \param parameters Parameter values to use for evaluation
     * \param variable Variable to return, if not specified then the variable specified in the constructor and stored in
     * var is returned
     *
     * If the same variables exist in parameters and in fixed_vals (set using set()), then the values from fixed_vals
     * are used.
     */
    template <class T>
    Mat<T> valueMat(const std::map<std::string,T> &parameters, const std::string &variable = std::string()) const;

    /*!
     * \brief Set variable in fixed_vals
     * \param name Name of variable to set
     * \param val Value to set
     */
    void set(const std::string &name, double val);

private:
    std::string eq;                          //!< Equation as string
    std::string var;                         //!< Standard return variable when value()/valueMat() is called
    std::map<std::string,double> fixed_vals; //!< Fixed values which were set using set()

    //! Stream operator for debug-output
    friend std::ostream & operator<<(std::ostream &os, const Equation& eq);
};

//! Stream operator for debug-output
std::ostream & operator<<(std::ostream &os, const Equation& eq);

#endif // EQUATION_H
