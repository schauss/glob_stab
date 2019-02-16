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

#include "parse_parameter.h"

#include <cstdio>
#include <boost/spirit/include/phoenix.hpp>

#include "parameters.h"
#include "iv.h"

using namespace boost::spirit;

ParseParameter::ParseParameter() : ParseParameter::base_type(start)
{
    using boost::phoenix::at_c;
    using boost::phoenix::assign;

    varName = raw[( qi::alpha | '_') >> *( qi::alnum | '_')];

    start =
        eps                       [at_c<3>(_val) = 0.0]
        >>  varName               [at_c<0>(_val) = _1]
        >> '(' >> double_         [at_c<1>(_val) = _1] [at_c<2>(_val) = _1]
        >> -( ',' >> double_      [at_c<2>(_val) = _1])
        >> -( ',' >> double_      [at_c<3>(_val) = _1])
        >> ')'
        >> *lit(';')
        ;
}

void ParseParameter::parse(const std::string &str)
{
    ParserVarParameter par;

    qi::rule<std::string::const_iterator> skip_ws_and_comments =
        qi::space | "#" >> *(qi::char_-qi::eol) >> qi::eol;

    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();

    bool r;
    do {
        r = phrase_parse( iter, end,*this,skip_ws_and_comments, par);
        if (r)
            Parameters::add(par.name,Iv(par.lower, par.upper),par.min_fraction);
    } while (r && (iter != end));

    if (!(r && (iter == end))) {
        std::string rest(iter, end);
        std::cout << "-------------------------\n";
        std::cout << "Parameter Parsing failed\n";
        std::cout << "stopped at: \"" << rest << "\"\n";
        std::cout << "-------------------------\n";
    }
}
