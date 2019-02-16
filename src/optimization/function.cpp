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

#include "function.h"

#include "config.h"
#include "value_set.h"

Function::Function(const Equation &eq, const Parameters &param)
    : param(param),
      evaluated(false),
      evaluated_from_equation(false),
      has_marginal(false)
{
    equation = boost::shared_ptr<const Equation>(new Equation(eq));
    std::cout << *equation << std::endl;
    std::cout << param.print() << std::endl;
}

Function::Function(const ValueSet &vs, const Parameters &param)
    : param(param),
      evaluated(false),
      evaluated_from_equation(false),
      has_marginal(false)
{
    value_set = boost::shared_ptr<const ValueSet>(new ValueSet(vs));
    std::cout << *value_set << std::endl;
    std::cout << param.print() << std::endl;
}

void Function::evaluateFromEquation()
{
    if (evaluated_from_equation)
        return;

    if (value_set.get()) {
        Vec<TaylorModel<Iv> > tmp = value_set->value(param.tmMap<Iv>());
        poly2d.clear();
        rest.clear();
        poly2d.push_back(Bernstein2d<Iv>(tmp));
        for (size_t i=0; i<tmp.size(); i++) {
            rest.push_back(tmp[i].restBound());
        }
    } else if (equation.get()) {
        Vec<TaylorModel<Iv> > coeffs = equation->value(param.tmMap<Iv>());
        poly.clear();
        rest.clear();
        for (size_t i=0; i<coeffs.size(); i++) {
            poly.push_back(Bernstein<Iv>(coeffs[i].polynomial_01()));
            rest.push_back(coeffs[i].restBound());
        }
    } else
        throw("Function::evaluateFromEquation called on empty function!");

    evaluated_from_equation = true;
}

void Function::evaluate()
{
    if (evaluated)
        return;

    if (poly.empty() && poly2d.empty()) {
        evaluated_from_equation = false;
        evaluateFromEquation();
    }

    if (!poly2d.empty())
        poly2d[0].evaluate();
    else {
        done = Vec<bool>(poly.size(),false);
        for (size_t i=0; i<poly.size(); i++) {
            poly[i].bound();
            if ( (poly[i].bound() > 0.0) || (poly[i].bound() <= 0.0) )
                done[i] = true;
        }
    }

    evaluated = true;
}

void Function::clear()
{
    poly.clear();
    poly2d.clear();
    evaluated = false;
}

Vec<Iv> Function::value() const
{
    if (!evaluated)
        throw("Function::value() - not evaluated!");

    if (value_set.get())
        throw("Function::value() not defined for value_set!");
    else
        return apply(poly,&Bernstein<Iv>::bound)+rest;
}

Vec<Vec<double> > Function::dpoly_dparam() const
{
    Vec<Vec<double> > ret;
    if (!poly2d.empty()) {
        ret = Vec<Vec<double> >(poly2d.size(),Vec<double>(param.size(),0.0));
        for (size_t i=0; i<poly2d.size(); i++) {
            ret[i] = poly2d[i].dx();
        }
    } else {
        ret = Vec<Vec<double> >(poly.size(),Vec<double>(param.size(),0.0));
        for (size_t i=0; i<poly.size(); i++) {
            if (done[i])
                ret[i] = Vec<double>(param.size(),0);
            else
                ret[i] = poly[i].dx();
        }
    }
    return ret;
}

Vec<Vec<double> > Function::drest_dparam() const
{
    Vec<Vec<double> > ret;
    if (value_set.get()) {
        bool use_rest = true;
        if ( (rest.size() == 2) && (rest[0].width()*rest[0].width() + rest[1].width()*rest[1].width() == 0.0) )
            use_rest = false;
        ret = Vec<Vec<double> >(1,Vec<double>(param.size(),0.0));
        for (size_t i=0; i<param.size(); i++) {
            std::vector<Parameters> tmp = param.splitEqual(i);
            Vec<TaylorModel<Iv> > tmp1 = value_set->value(tmp[0].tmMap<Iv>());

            double tmp2;
            if (use_rest)
                tmp2 = tmp1[0].restBound().width()*tmp1[0].restBound().width()+tmp1[1].restBound().width()*tmp1[1].restBound().width();
            else
                tmp2 = tmp1[0].bound().width()*tmp1[0].bound().width()+tmp1[1].bound().width()*tmp1[1].bound().width();

            ret[0][i] = tmp2;
        }
    } else {
        Vec<Vec<double> > dx;
        for (size_t i=0; i<param.size(); i++) {
            std::vector<Parameters> tmp = param.splitEqual(i);
            Vec<TaylorModel<Iv> > tmp0 = equation->value(tmp[0].tmMap<Iv>());
            Vec<TaylorModel<Iv> > tmp1 = equation->value(tmp[1].tmMap<Iv>());
            Vec<double> tmp2;
            for (size_t j=0;j<tmp1.size();j++) {
                if (done[j])
                    tmp2.push_back(-INFINITY);
                else {
                    //tmp2.push_back(tmp1[j].restBound().width());
                    double v0 = tmp0[j].bound().boundViolation(0.0);
                    double v1 = tmp1[j].bound().boundViolation(0.0);
                    tmp2.push_back(v0>v1?v0:v1);
                    //tmp2.push_back((Bernstein<Iv>(tmp1[j].polynomial_01()).bound()+tmp1[j].restBound()).boundViolation(0.0));
                }
            }
            dx.push_back(tmp2);
        }

        ret = Vec<Vec<double> >(dx[0].size(),Vec<double>(param.size(),0.0));
        for (size_t i=0; i<ret.size(); i++) {
            for (size_t j=0; j<ret[0].size(); j++) {
                ret[i][j] = dx[j][i];
            }
        }
    }
    return ret;
}

bool Function::includesZero() const
{
    if (!evaluated)
        throw("Function::includesZero() - not evaluated!");

    if (!poly2d.empty()) {
        return poly2d[0].restIncluded();
    } else {
        for (size_t i=0; i<poly.size(); i++) {
            if ( poly[i].innerBound().includes(-rest[i]) )
                return true;
        }
        return false;
    }
}

bool Function::partiallyIncludesZero1() const
{
    if (!evaluated)
        throw("Function::partiallyIncludesZero1() - not evaluated!");

    if (!poly2d.empty()) {
        return ( poly2d[0].restPartiallyIncluded() );
    } else {
        for (size_t i=0; i<poly.size(); i++) {
            //if (!poly[i].innerBound().includes(0.0) && poly[i].bound().includes(0.0) )
            if (!poly[i].innerBound().includes(-rest[i]) && poly[i].bound().includes(-rest[i]) )
                return true;
        }
        return false;
    }
}

bool Function::partiallyIncludesZero2() const
{
    if (!evaluated)
        throw("Function::partiallyIncludesZero2() - not evaluated!");

    if (!poly2d.empty()) {
        return ( poly2d[0].restPartiallyExcluded() || poly2d[0].zeroIncluded() );
    } else {
        for (size_t i=0; i<poly.size(); i++) {
            if (!poly[i].innerBound().includes(-rest[i]) && (poly[i].bound()+rest[i]).includes(0.0) )
                return true;
        }
        return false;
    }
}

bool Function::excludesZero() const
{
    if (!evaluated)
        throw("Function::excludesZero() - not evaluated!");

    if (!poly2d.empty()) {
        return poly2d[0].restExcluded();
    } else {
        if ( (value()<=0.0).any() || (value()>0.0).all())
            return true;
        else
            return false;
    }
}

bool Function::boundaryIsZero() const
{
    if (!evaluated)
        throw("Function::boundaryIsZero() - not evaluated!");

    if (!poly2d.empty())
        return poly2d[0].boundaryIsZero();
    else {
        Vec<Iv> val = value();
        for (size_t i=0; i<poly.size(); i++) {
            if ( (poly[i].innerBound().lower() == 0.0) || (poly[i].innerBound().upper() == 0.0) ||
                 (poly[i].bound().lower() == 0.0) || (poly[i].bound().upper() == 0.0) ||
                 (val[i].lower() == 0.0) || (val[i].upper() == 0.0) )
                return true;
        }
        return false;
    }
}

bool Function::restIsFinite() const
{
    if (!evaluated)
        throw("Function::restIsFinite() - not evaluated!");

    for (size_t i=0; i<rest.size(); i++) {
        if (!rest[i].isFinite())
            return false;
    }
    return true;
}

std::vector<boost::shared_ptr<Function> > Function::splitAt(size_t direction, double mid_point) const
{
    std::vector<boost::shared_ptr<Function> > ret;
    ret.push_back(boost::shared_ptr<Function>(new Function(*this)));
    ret.push_back(boost::shared_ptr<Function>(new Function(*this)));

    // generate the parameter sets
    Vec<Parameters> param_new = param.splitAtRelative(direction,mid_point);

    for(size_t i=0; i<ret.size(); i++) {
        ret[i]->evaluated = false;
        ret[i]->evaluated_from_equation = false;
        ret[i]->param = param_new[i];
    }

    for(size_t i=0; i<poly.size(); i++) {
        Vec<Bernstein<Iv> > poly_new = poly[i].splitAt(direction, mid_point);
        for (size_t j=0; j<ret.size(); j++)
            ret[j]->poly[i] = poly_new[j];
    }

    for(size_t i=0; i<poly2d.size(); i++) {
        Vec<Bernstein2d<Iv> > poly2d_new = poly2d[i].splitAt(direction, mid_point);
        for (size_t j=0; j<ret.size(); j++)
            ret[j]->poly2d[i] = poly2d_new[j];
    }

    return ret;
}

void Function::enlarge()
{
    clear();
    param.scale(1.0+1e-10);
}

bool Function::hasMarginal()
{
    return has_marginal;
}

void Function::setMarginal(bool has_marginal)
{
    this->has_marginal = has_marginal;
}

std::ostream & operator<<(std::ostream &os, Function& fcn)
{
    //fcn.evaluate();

    std::ostringstream buf;

    buf << fcn.param.print();
    buf << "poly-size: " << fcn.poly.size() << ", poly2d-size: " << fcn.poly2d.size() << ", rest-size: " << fcn.rest.size() << std::endl;
    if (!fcn.poly.empty()) {
        buf << "Value:           " << fcn.value() << std::endl;
        buf << "Poly-Value:      " << apply(fcn.poly,&Bernstein<Iv>::bound) << std::endl;
        buf << "Poly-Value-inner:" << apply(fcn.poly,&Bernstein<Iv>::innerBound) << std::endl;
        buf << "Rest:            " << fcn.rest << std::endl;
        if (CONFIG_VAR(verbosity,int) > 1) {
            for(size_t i=0; i<fcn.poly.size(); i++) {
                buf << "poly[" << i << "]: " << fcn.poly[i].coeffs() << std::endl;
            }
        }
    } else if (!fcn.poly2d.empty()) {
        if (CONFIG_VAR(verbosity,int) > 1) {
            for(size_t i=0; i<fcn.poly2d.size(); i++) {
                buf << "poly2d[" << i << "]: " << std::endl << fcn.poly2d[i] << std::endl;
            }
        }
    }

    if (fcn.doValueSet()) {
        buf << "VS-Boundary:     " << fcn.value_set->doBoundary() << std::endl;
        buf << "VS-Inverse:      " << fcn.value_set->doInverse() << std::endl;
    }

    return os << buf.str();
}
