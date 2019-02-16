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

#include "parse_identifier.h"
#include <cstdio>

#include <boost/spirit/include/phoenix.hpp>

using namespace boost::spirit;

ParseIdentifier::ParseIdentifier() : ParseIdentifier::base_type(start)
{
    using boost::phoenix::at_c;
    using boost::phoenix::assign;

    varName = raw[( qi::alpha | '_') >> *( qi::alnum | '_')];

    start =
        eps [at_c<1>(_val) = -1] [at_c<2>(_val) = -1] [at_c<3>(_val) = -1] [at_c<4>(_val) = -1]
        >>  varName             [at_c<0>(_val) = _1]
        >> -( '[' >> int_       [at_c<1>(_val) = _1]
        >> -( ':' >> int_       [at_c<2>(_val) = _1] )
        >> -( ',' >> int_       [at_c<3>(_val) = _1]
        >> -( ':' >> int_       [at_c<4>(_val) = _1] ) )
        >> ']' )
        ;
}

ParserVarId ParseIdentifier::parse(const std::string &str) const
{
    ParserVarId id;
    boost::spirit::ascii::space_type space;

    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bool r = phrase_parse( iter, end,*this,space, id);

    if (!(r && (iter == end))) {
        std::string rest(iter, end);
        std::cout << "-------------------------\n";
        std::cout << "ID Parsing failed\n";
        std::cout << "stopped at: \"" << rest << "\"\n";
        std::cout << "-------------------------\n";
    }

    return id;
}
