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

#include "equation.h"

#include "parse_equation.h"
#include "parameters.h"
#include "taylor_model.h"

Equation::Equation(const std::string &equation, const std::string &variable)
    : eq(equation), var(variable)
{

}

template <class T>
Vec<T> Equation::value(const std::map<std::string,T> &parameters, const std::string &variable) const
{
    ParseEquation<T> parser;
    parser.set(parameters);

    for (typename std::map<std::string,double>::const_iterator it=fixed_vals.begin(); it != fixed_vals.end(); it++)
        parser.set(it->first,T(it->second));

    parser.parse(this->eq);

    if (variable.empty())
        return parser.getVec(this->var);
    else
        return parser.getVec(variable);
}

template Vec<double>               Equation::value(const std::map<std::string,double>               &parameters, const std::string &variable) const;
template Vec<TaylorModel<double> > Equation::value(const std::map<std::string,TaylorModel<double> > &parameters, const std::string &variable) const;
template Vec<TaylorModel<Iv> >     Equation::value(const std::map<std::string,TaylorModel<Iv> >     &parameters, const std::string &variable) const;
template Vec<std::complex<double> >               Equation::value(const std::map<std::string,std::complex<double> >               &parameters, const std::string &variable) const;
template Vec<std::complex<TaylorModel<double> > > Equation::value(const std::map<std::string,std::complex<TaylorModel<double> > > &parameters, const std::string &variable) const;
template Vec<std::complex<TaylorModel<Iv> > >     Equation::value(const std::map<std::string,std::complex<TaylorModel<Iv> > >     &parameters, const std::string &variable) const;

template <class T>
Mat<T> Equation::valueMat(const std::map<std::string,T> &parameters, const std::string &variable) const
{
    ParseEquation<T> parser;
    parser.set(parameters);

    for (typename std::map<std::string,double>::const_iterator it=fixed_vals.begin(); it != fixed_vals.end(); it++)
        parser.set(it->first,T(it->second));

    parser.parse(this->eq);

    if (variable.empty())
        return parser.getMat(this->var);
    else
        return parser.getMat(variable);
}

template Mat<double>               Equation::valueMat(const std::map<std::string,double>               &parameters, const std::string &variable) const;
template Mat<TaylorModel<double> > Equation::valueMat(const std::map<std::string,TaylorModel<double> > &parameters, const std::string &variable) const;
template Mat<TaylorModel<Iv> >     Equation::valueMat(const std::map<std::string,TaylorModel<Iv> >     &parameters, const std::string &variable) const;
template Mat<std::complex<double> >               Equation::valueMat(const std::map<std::string,std::complex<double> >               &parameters, const std::string &variable) const;
template Mat<std::complex<TaylorModel<double> > > Equation::valueMat(const std::map<std::string,std::complex<TaylorModel<double> > > &parameters, const std::string &variable) const;
template Mat<std::complex<TaylorModel<Iv> > >     Equation::valueMat(const std::map<std::string,std::complex<TaylorModel<Iv> > >     &parameters, const std::string &variable) const;

void Equation::set(const std::string &name, double val)
{
    fixed_vals[name] = val;
}

std::ostream & operator<<(std::ostream &os, const Equation& eq)
{
    std::ostringstream buf;

    buf << "Equation:\n" << eq.eq << std::endl;
    buf << "Variable:\n" << eq.var << std::endl;

    return (os << buf.str());
}
