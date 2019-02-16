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

#ifndef PARSEPARAMETER_H
#define PARSEPARAMETER_H

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <string>

using boost::spirit::qi::symbols;

//! The ParserVarParameter struct describes a parameter
struct ParserVarParameter
{
    std::string name;    //!< Name of the parameter
    double lower;        //!< Lower bound of the parameter
    double upper;        //!< Upper bound of the parameter
    double min_fraction; //!< Relative resolution of boundary mapping in direction of this parameter
};

/*!
 * \brief The class ParseParameter is used to parse a uncertain parameter.
 *
 * This can be something like:
 * * var(1.0) => Iv(1.0,1.0)
 * * var(1.0,2.0) => Iv(1.0,2.0)
 * * var(1.0,2.0,0.1) )> Iv(1.0,2.0), resolution of boundary mapping below 0.1*(2.0-1.0) in this dimension
 *
 * The variable name var can be something like: a, _a, b0, b1, c123, A, _Bbsdf234
 */
class ParseParameter : public boost::spirit::qi::grammar<std::string::const_iterator, ParserVarParameter(), boost::spirit::qi::rule<std::string::const_iterator> >
{
public:
    //! Instantiate ParseIdentifier and define grammar
    ParseParameter();
    //! Parse string and add parameters to Parameters
    void parse(const std::string &str);

private:
    boost::spirit::qi::rule<std::string::const_iterator, std::string()> varName;
    boost::spirit::qi::rule<std::string::const_iterator, ParserVarParameter(), boost::spirit::qi::rule<std::string::const_iterator> > start;
};

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

BOOST_FUSION_ADAPT_STRUCT(
    ParserVarParameter,
    (std::string, name)
    (double, lower)
    (double, upper)
    (double, min_fraction)
)

#endif // PARSEPARAMETER_H
