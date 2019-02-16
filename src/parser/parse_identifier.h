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

#ifndef PARSEIDENTIFIER_H
#define PARSEIDENTIFIER_H

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <string>

using boost::spirit::qi::symbols;

//! The ParserVarId struct describes a Variable
struct ParserVarId
{
    std::string name; //!< Name of the variable
    int start_row;    //!< Start row (-1 for scalar)
    int end_row;      //!< End row (last element + 1, -1 for scalar or element of vector / row of matrix, i.e., only start row is given)
    int start_col;    //!< Start column (-1 for scalar or vector)
    int end_col;      //!< End column (last element + 1, -1 for scalar or vector)
};

/*!
 * \brief The class ParseIdentifier is used to parse a variable identifier within an equation.
 *
 * This can be something like:
 * * a, _a, b0, b1, c123 (scalar)
 * * a[0], a[1], a[2] (vector)
 * * a[0:3] (elements 0,1,2 of vector)
 * * B[0,0], B[1,0], B[1,1] (matrix)
 * * B[1:3,2:4] (elements B[1,2],B[1,3],B[2,2],B[2,3] of matrix B)
 */
class ParseIdentifier : public boost::spirit::qi::grammar<std::string::const_iterator, ParserVarId()>
{
public:
    //! Instantiate ParseIdentifier and define grammar
    ParseIdentifier();
    //! Parse string returning ParserVarId
    ParserVarId parse(const std::string &str) const;

private:
    boost::spirit::qi::rule<std::string::const_iterator, std::string()> varName;
    boost::spirit::qi::rule<std::string::const_iterator, ParserVarId()> start;
};

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

BOOST_FUSION_ADAPT_STRUCT(
    ParserVarId,
    (std::string, name)
    (int, start_row)
    (int, end_row)
    (int, start_col)
    (int, end_col)
)

#endif // PARSEIDENTIFIER_H
