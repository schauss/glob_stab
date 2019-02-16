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

#include "parse_equation.h"
#include <cstdio>

#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/phoenix.hpp>
#include <stdexcept>

using namespace boost::spirit;

/*!
 * See the documentation of boost::spirit and especially qi and phoenix for details.
 */
template <typename T>
ParseEquation<T>::ParseEquation() : ParseEquation<T>::base_type(expression)
{
    expression =
        ( (identifier >> '=' >> sum)      [_val = boost::phoenix::bind(&ParseEquation<T>::set, this, _1, _2)]
          | sum                           [_val = _1] )
        >> *lit(';')
        ;

    identifier =
        raw[( qi::alpha | '_') >> *( qi::alnum | '_') >> -( '[' >> uint_ >> -( ':' >> uint_ ) >> -( ',' >> uint_ >> -( ':' >> uint_ ) ) >> ']')]
        ;
    sum =
        term                            [_val  = _1]
        >> *(   ('+' >> term            [_val += _1])
            |   ('-' >> term            [_val -= _1])
            )
        ;

    term =
        factor                          [_val  = _1]
        >> *(   ('*' >> factor          [_val *= _1])
            |   ('/' >> factor          [_val /= _1])
            )
        ;

    factor =
            ( lit("sin") >> '(' >> sum  [_val = boost::phoenix::bind(&ParseEquation<T>::intrinsic, this, _1,1)] >> ')' )
        |   ( lit("cos") >> '(' >> sum  [_val = boost::phoenix::bind(&ParseEquation<T>::intrinsic, this, _1,2)] >> ')' )
        |   ( lit("exp") >> '(' >> sum  [_val = boost::phoenix::bind(&ParseEquation<T>::intrinsic, this, _1,3)] >> ')' )
        |   identifier                  [_val = boost::phoenix::bind(&ParseEquation<T>::get, this, _1)]
        |   double_                     [_val = _1]
        |   '(' >> sum                  [_val = _1] >> ')'
        |   ('-' >> factor              [_val = -_1])
        |   ('+' >> factor              [_val =_1])
        ;
}

template <typename T>
T ParseEquation<T>::set (const std::string& name, const T &value)
{
    ParserVarId id = idParser.parse(name);
    return setVar(id,value);
}

template <typename T>
void ParseEquation<T>::set (const std::map<std::string,T> &value)
{
    for (typename std::map<std::string,T>::const_iterator it=value.begin(); it!= value.end(); it++) {
        ParserVarId id = idParser.parse(it->first);
        setVar(id,it->second);
    }
}

template <typename T>
T ParseEquation<T>::get (const std::string& name) const
{
    ParserVarId id = idParser.parse(name);
    Mat<T> v = getVar(id);
    if ( (v.size() == 1) && (v[0].size() == 1) )
        return v[0][0];
    else
        throw("ParseEquation<T>::get - variable is no scalar!");
}

template <typename T>
Vec<T> ParseEquation<T>::getVec (const std::string& name) const
{
    ParserVarId id = idParser.parse(name);
    Mat<T> v = getVar(id);
    if (v.size() == 1)
        return v[0];
    else {
        Vec<T> ret;
        for (size_t i = 0; i<v.size(); i++) {
            if (v[i].size() != 1) {
                std::cout << "Error in ParseEquation<T>::getVec(" << name << "): row[" << i << "].size() = " << v[i].size() << std::endl;
                throw("ParseEquation<T>::getVec - variable is no vector!");
            } else {
                ret.push_back(v[i][0]);
            }
        }
        return ret;
    }
}

template <typename T>
Mat<T> ParseEquation<T>::getMat (const std::string& name) const
{
    ParserVarId id = idParser.parse(name);
    return getVar(id);
}

template <typename T>
T ParseEquation<T>::intrinsic (const T &value, size_t type) const
{
    if (type==1)
        return sin(value);
    else if (type==2)
        return cos(value);
    else if (type==3)
        return exp(value);
    else
        throw("ParseEquation<T>::intrinsic - unknown function!");
}

/*
template <typename T>
T ParseEquation<T>::cast (const double &value) const
{
    return Vec<T>(1,(T)value);
}
*/

template <typename T>
T ParseEquation<T>::setVar(const ParserVarId id, T value)
{
    if (id.start_row < 0) {
        var[id.name] = Mat<T>();
        var[id.name].push_back(Vec<T>());
        var[id.name][0].push_back(Vec<T>(1,value));
        return value;
    }

    typename std::map<std::string,Mat<T> >::iterator iter = var.find(id.name);
    if (iter == var.end()) {
        var[id.name] = Mat<T>();
        iter = var.find(id.name);
    }

    int start_row = id.start_row>=0?id.start_row:0;
    int end_row = id.end_row>0?id.end_row:start_row+1;
    int start_col = id.start_col>=0?id.start_col:0;
    int end_col = id.end_col>0?id.end_col:start_col+1;

    while ((int)iter->second.size() < end_row) {
        iter->second.push_back(Vec<T>());
    }

    for (size_t i = (size_t)start_row; i<(size_t)end_row; i++) {
        if (iter->second.size() <= i)
            iter->second.push_back(Vec<T>());
        for (size_t j = (size_t)start_col; j<(size_t)end_col; j++) {
            if (iter->second[i].size() <= j)
                iter->second[i].push_back(T());
            iter->second[i][j] = value;
        }
    }

    return value;
}

template <typename T>
Mat<T> ParseEquation<T>::getVar(const ParserVarId id) const
{
    typename std::map<std::string,Mat<T> >::const_iterator iter = var.find(id.name);
    if (iter == var.end()) {
        std::cerr << __PRETTY_FUNCTION__ << "Could not find variable " << id.name << "!\n";
        throw("Parsing error!");
    }

    if (id.start_row < 0) {
        return (iter->second);
    }

    //std::cout << "Called getVar: " << id.name << "[" << id.start_row << ":" << id.end_row << ", " << id.start_col << ":" << id.end_col << "]" << std::endl;

    int start_row = id.start_row>=0?id.start_row:0;
    int end_row = id.end_row>0?id.end_row:(id.start_row>=0?start_row+1:iter->second.size());
    int start_col = id.start_col>=0?id.start_col:0;
    int end_col = id.end_col>0?id.end_col:(id.start_col>=0?start_col+1:iter->second[start_row].size());

    Mat<T> ret;
    for (size_t i = (size_t)start_row; i<(size_t)end_row; i++) {
        ret.push_back(Vec<T>());
        for (size_t j = (size_t)start_col; j<(size_t)end_col; j++) {
            ret.back().push_back(iter->second[i][j]);
        }
    }

    return ret;
}

template <typename T>
void ParseEquation<T>::clear(std::string name){
    if (name.empty())
        var.clear();
    else {
        typename std::map<std::string,Mat<T> >::iterator iter = var.find(name);
        if (iter == var.end()) {
            std::cerr << __PRETTY_FUNCTION__ <<  " - Could not find variable " << name << "!\n";
        }
        else
            var.erase(iter);
    }

}

template <typename T>
T ParseEquation<T>::parse(const std::string &str)
{
    qi::rule<std::string::const_iterator> skip_ws_and_comments =
        qi::space | "#" >> *(qi::char_-qi::eol) >> qi::eol;

    if (str.empty()) {
        throw(std::runtime_error(std::string(__PRETTY_FUNCTION__)+": String empty!\n"));
    }

    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    T result;
    bool r;
    do {
        r = phrase_parse( iter, end, *this, skip_ws_and_comments, result);
    } while (r && (iter != end));

    if (!(r && (iter == end))) {
        std::string rest(iter, end);
        std::cout << "-------------------------\n";
        std::cout << "Parsing failed\n";
        std::cout << "stopped at: \"" << rest << "\"\n";
        std::cout << "-------------------------\n";
        throw(std::runtime_error(std::string(__PRETTY_FUNCTION__)+": Parsing failed!\n"));
    }

    return result;
}

#include "taylor_model.h"
template class ParseEquation<int>;
template class ParseEquation<double>;
template class ParseEquation<TaylorModel<double> >;
template class ParseEquation<TaylorModel<Iv> >;
template class ParseEquation<std::complex<double> >;
template class ParseEquation<std::complex<TaylorModel<double> > >;
template class ParseEquation<std::complex<TaylorModel<Iv> > >;
