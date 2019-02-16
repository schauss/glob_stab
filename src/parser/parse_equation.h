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

#ifndef PARSEEQUATION_H
#define PARSEEQUATION_H

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <string>
#include <map>
#include "vec.hpp"
#include "mat.hpp"
#include "parse_identifier.h"

using boost::spirit::qi::symbols;

/*!
 * \brief The class ParseEquation is used to evaluate a mathematical expression for a certain type.
 *
 * The class is used within this program to evaluate the characteristic function of a dynamic system using Taylor
 * Models. The implementation is based on boost::spirit::qi which allows to simply define a grammar and parse
 * expressions.
 *
 * \todo As standard data type Mat<T> is used. This rsults in quite some overhead as we are often only working with
 * scalars or at the very most vectors => optimization possible.
 */
template <typename T>
class ParseEquation : public boost::spirit::qi::grammar<std::string::const_iterator, T(), boost::spirit::qi::rule<std::string::const_iterator> >
{
public:
    /*!
     * \brief Within the constructor the grammar is defined.
     */
    ParseEquation();

    //! Parse the string str returning the result (which is the value of the last expression)
    T parse(const std::string &str);

    /*!
     * \brief Set a variable to a value
     * \param name Variable name
     * \param value New variable value
     * \return Returns the value to which the variable was set
     */
    T      set  (const std::string& name, const T &value);
    /*!
     * \brief Set several variables to given values
     * \param value map with variable name as key and variable value as value
     */
    void   set  (const std::map<std::string,T> &value);

    /*!
     * \brief Return a scalar variable with given name
     * \param name Variable name
     */
    T      get     (const std::string& name) const;
    /*!
     * \brief Return a vector variable with given name
     * \param name Variable name
     */
    Vec<T> getVec  (const std::string& name) const;
    /*!
     * \brief Return a matrix variable with given name
     * \param name Variable name
     */
    Mat<T> getMat  (const std::string& name) const;

    //! Return all variables as map
    std::map<std::string,Mat<T> > variables() const {return var;}

    /*!
     * \brief Remove existing variables
     * \param name Variable name, if empty then all variables are cleared
     */
    void clear(std::string name = std::string());

private:
    /*!
     * \brief Evaluate intrinsic function
     * \param value Value of argument
     * \param type 1: sin, 2: cos, 3: exp
     * \return Result of evaluating intrinsic function with value
     */
    T intrinsic (const T &value, size_t type) const;
    //T cast (const double &value) const;

    boost::spirit::qi::rule<std::string::const_iterator, T(), boost::spirit::qi::rule<std::string::const_iterator> > expression, term, factor, sum;
    boost::spirit::qi::rule<std::string::const_iterator, std::string()> identifier;

    //! All variables are stored in this map. This includes parameters as well as all results while evaluating the equation.
    std::map<std::string,Mat<T> > var;

    /*!
     * \brief Set Variable with given id to given value
     * \param id Variable id
     * \param value New Value
     * \return New value
     *
     * This can also be used to set a whole matrix or a rectangular submatrix to the value
     */
    T setVar(const ParserVarId id, T value);

    /*!
     * \brief Return variable with given id
     * \param id Variable id
     */
    Mat<T> getVar(const ParserVarId id) const;
    //typename std::map<std::string,Mat<T> >::iterator getIter (const std::string &str);

    ParseIdentifier idParser;
};

#endif // PARSEEQUATION_H
