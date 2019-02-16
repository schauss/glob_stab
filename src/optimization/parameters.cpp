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

#include "parameters.h"
#include "parse_parameter.h"

#include <cstdio>

Parameters          Parameters::original            = Parameters(true);
Vec<std::string>    Parameters::var_name            = Vec<std::string>();
Vec<std::string>    Parameters::var_fixed_name      = Vec<std::string>();
Vec<double>         Parameters::min_fraction        = Vec<double>();
bool                Parameters::fraction_set        = false;
std::map<std::string,size_t> Parameters::var_id     = std::map<std::string,size_t>();

Parameters::Parameters() :
    val(original.val),
    val_fixed(original.val_fixed),
    start_idx(original.start_idx),
    end_idx(original.end_idx)
{
    if (original.val.empty())
        throw("Parameters::Parameters(): Initial parameters empty!");
}

/*!
 * \bug Beware, as this uses idxToD it may underapproximate the result due to
 * rounding errors!
 */
Parameters::Parameters(const Vec<size_t> &start, const Vec<size_t> &end) :
    val_fixed(original.val_fixed),
    start_idx(start),
    end_idx(end)
{
    if ( (start.size() != original.start_idx.size()) ||
         (end.size() != original.end_idx.size()) )
        throw("Parameters::Parameters(start, end) - index-size missmatch!");

    for (size_t i=0; i<start.size(); i++) {
        if ( (start[i] < original.start_idx[i]) || (end[i] > original.end_idx[i]) )
            throw("Parameters::Parameters(start, end) - index out of bounds!");

        val.push_back( Iv(original.idxToD(start[i],i,false),original.idxToD(end[i],i,false)) );
    }
}

/*!
 * Necessary as the default constructor Paramters() actually initializes several
 * member variables to the same value as member variables of
 * Parameters::original. So we can obviously not call that constructor to create
 * Paramters::original in the first place.
 */
Parameters::Parameters(bool stat)
{
}

/*!
 * Values stored in Parameters::val as well as Parameters::val_fixed are returned.
 */
Vec<Iv> Parameters::value() const {
    Vec<Iv> ret = val;
    ret.push_back(val_fixed);
    return ret;
}

/*!
 * Values stored in Parameters::val as well as Parameters::val_fixed are returned.
 */
template <class T_coeff>
Vec<TaylorModel<T_coeff> > Parameters::tmValue() const {
    Vec<TaylorModel<T_coeff> > ret;
    for (size_t i=0; i<val.size(); i++) {
        ret.push_back(TaylorModel<T_coeff>(var_name[i],val[i]));
    }
    for (size_t i=0; i<val_fixed.size(); i++) {
        ret.push_back(TaylorModel<T_coeff>(val_fixed[i]));
    }
    return ret;
}

template Vec<TaylorModel<double> > Parameters::tmValue<double>() const;
template Vec<TaylorModel<Iv> > Parameters::tmValue<Iv>() const;

/*!
 * Values stored in Parameters::val as well as Parameters::val_fixed are returned.
 */
template <class T_coeff>
Vec<std::complex<TaylorModel<T_coeff> > > Parameters::tmValueComplex() const {
    Vec<std::complex<TaylorModel<T_coeff> > > ret;
    for (size_t i=0; i<val.size(); i++) {
        ret.push_back(std::complex<TaylorModel<T_coeff> >(TaylorModel<T_coeff>(var_name[i],val[i]),false));
    }
    for (size_t i=0; i<val_fixed.size(); i++) {
        ret.push_back(std::complex<TaylorModel<T_coeff> >(TaylorModel<T_coeff>(val_fixed[i]),false));
    }
    return ret;
}

template Vec<std::complex<TaylorModel<double> > > Parameters::tmValueComplex<double>() const;
template Vec<std::complex<TaylorModel<Iv> > > Parameters::tmValueComplex<Iv>() const;

/*!
 * Values stored in Parameters::val as well as Parameters::val_fixed are returned.
 */
std::map<std::string,Iv> Parameters::map() const {
    std::map<std::string,Iv> ret;
    for (size_t i=0; i<val.size(); i++) {
        ret[var_name[i]] = val[i];
    }
    for (size_t i=0; i<val_fixed.size(); i++) {
        ret[var_fixed_name[i]] = val_fixed[i];
    }
    return ret;
}

/*!
 * Values stored in Parameters::val as well as Parameters::val_fixed are returned.
 */
template <class T_coeff>
std::map<std::string,TaylorModel<T_coeff> > Parameters::tmMap() const
{
    std::map<std::string,TaylorModel<T_coeff> > ret;
    for (size_t i=0; i<val.size(); i++) {
        ret[var_name[i]] = TaylorModel<T_coeff>(var_name[i],val[i]);
    }
    for (size_t i=0; i<val_fixed.size(); i++) {
        ret[var_fixed_name[i]] = TaylorModel<T_coeff>(val_fixed[i]);
    }
    return ret;
}

template std::map<std::string,TaylorModel<double> > Parameters::tmMap<double>() const;
template std::map<std::string,TaylorModel<Iv> > Parameters::tmMap<Iv>() const;

/*!
 * Values stored in Parameters::val as well as Parameters::val_fixed are returned.
 */
template <class T_coeff>
std::map<std::string,std::complex<TaylorModel<T_coeff> > > Parameters::tmMapComplex() const
{
    std::map<std::string,std::complex<TaylorModel<T_coeff> > > ret;
    for (size_t i=0; i<val.size(); i++) {
        ret[var_name[i]] = std::complex<TaylorModel<T_coeff> >(TaylorModel<T_coeff>(var_name[i],val[i]),false);
    }
    for (size_t i=0; i<val_fixed.size(); i++) {
        ret[var_fixed_name[i]] = std::complex<TaylorModel<T_coeff> >(TaylorModel<T_coeff>(val_fixed[i]),false);
    }
    return ret;
}

template std::map<std::string,std::complex<TaylorModel<double> > > Parameters::tmMapComplex<double>() const;
template std::map<std::string,std::complex<TaylorModel<Iv> > > Parameters::tmMapComplex<Iv>() const;

size_t Parameters::size() const
{
    return original.val.size();
}

/*!
 * First the parameter set is copied twice and the interval with given \a
 * direction is split at \a mid_point. Then, the split point is actually
 * adjusted to lie exactly on the boarder of two steps at desired minimum
 * resolution if a desired minimum resolution is given for this dimension and
 * this resolution has not been reached yet.
 *
 * When splitting into two equal parts this can be indicated by setting \a
 * middle to \em true which allows an easier calculation of the exact split-
 * point.
 * \todo mid_point is not a good name, better split_point?
 * \todo middle is not a good name for the parameter, moreover we should be able
 * to get rid of this parameter and do this automatically! Or add a function
 * splitEqual().
 */
Vec<Parameters> Parameters::splitAt(size_t direction, double mid_point, bool middle) const
{
    if ( (mid_point < val[direction].lower()) || (mid_point > val[direction].upper()) ) {
        std::cout << "lower, midpoint, upper: " << val[direction].lower() << ", " << mid_point << ", " << val[direction].upper() << std::endl;
        std::cout << print() << std::endl;
        throw("Parameters::splitAt - out of range!\n");
    }

    Vec<Iv> param_split = val[direction].splitAt(mid_point);
    Vec<Parameters> ret;
    for (size_t i = 0; i < param_split.size(); i++) {
        Parameters tmp = *this;
        tmp.val[direction] = param_split[i];
        ret.push_back(tmp);
    }

    int idx_range = end_idx[direction]-start_idx[direction];
    if (idx_range > 1) {
        if (middle && ( idx_range%2 == 0) ) {
            ret[0].end_idx[direction] =   start_idx[direction]+idx_range/2;
            ret[1].start_idx[direction] = start_idx[direction]+idx_range/2;
        } else {
            ret[0].end_idx[direction] = original.dToIdx(mid_point,direction,true);
            ret[1].start_idx[direction] = original.dToIdx(mid_point,direction,false);
        }
    }

    //std::cout << "idx_mid = " << idx_mid << std::endl;

    return ret;
}

/*!
 * Calculates the absolute split point and then calls
 * splitAt(size_t direction, double mid_point, bool middle) const where \a
 * middle is set to \em true if \a mid_point is \em 0.5.
 */
Vec<Parameters> Parameters::splitAtRelative(size_t direction, double mid_point) const
{
    double abs_mid = val[direction].lower() + mid_point * val[direction].width();
    return splitAt(direction,abs_mid,(mid_point==0.5));
}

/*!
 * Calculates the absolute split points and then calls
 * splitAt(size_t direction, double mid_point, bool middle) const as often as
 * necessary always adding the left parameter set and adding the right set for
 * the last split. If \a num_parts is \em 2 then \a middle is set to \em true.
 */
Vec<Parameters> Parameters::splitEqual(size_t direction, size_t num_parts) const
{
    Vec<Parameters> ret;

    double step = val[direction].width()/(double)num_parts;
    double split_point = val[direction].lower() + step;

    Parameters tmp = *this;
    for (size_t i = 0; i < num_parts-1; i++) {
        Vec<Parameters> split = tmp.splitAt(direction, split_point,(num_parts==2));
        ret.push_back(split[0]);
        tmp = split[1];
        split_point += step;
    }
    ret.push_back(tmp);

    return ret;
}

Parameters Parameters::intersect(const Parameters &p) const
{
    Parameters ret(*this);
    for (size_t i=0; i<size(); i++) {
        ret.val[i] = ret.val[i].intersect(p.val[i]);
        if (p.start_idx[i] > start_idx[i])
            ret.start_idx[i] = p.start_idx[i];
        if (p.end_idx[i] < end_idx[i])
            ret.end_idx[i] = p.end_idx[i];

        if ( val[i].isNaN() || (ret.start_idx[i] >= ret.end_idx[i]) )
            throw("Parameters::intersect - no intersection!");
    }
    return ret;
}

double Parameters::fraction(size_t i) const
{
    return (val[i].width()/original.val[i].width());
}

double Parameters::minFraction(size_t i) const
{
    return (min_fraction[i]);
}

bool Parameters::minFractionReached(int i) const
{
    if (!fraction_set)
        return true;

    if (i<0) {
        for (size_t j=0; j<size(); j++) {
            if (!minFractionReached(j))
                return false;
        }
        return true;
    }
    else {
        //return ( (minFraction((size_t)i)==0) || (fraction((size_t)i) < minFraction((size_t)i)) );
        return ( (minFraction((size_t)i)==0) || (end_idx[(size_t)i]-start_idx[(size_t)i] == 1) );
    }
}

size_t  Parameters::maxRelFractionDir() const
{
    int dir = -1;
    double rel_frac = 0.0;
    for (size_t i=0; i<size(); i++) {
        if (minFraction(i) == 0.0)
            continue;

        double tmp = fraction(i) / minFraction(i);
        if (tmp > rel_frac) {
            dir = i;
            rel_frac = tmp;
        }
    }

    if (dir == -1)
        throw("Parameters::maxRelFractionDir called for parameters where no fraction was set!");

    return (size_t)dir;
}

bool Parameters::isOriginal(size_t i) const
{
    return (val[i] == original.val[i]);
}

void Parameters::correct()
{
    for (size_t i = 0; i<size(); i++)
        val[i] = val[i].intersect(original.val[i]);
}

/*!
 * Does not affect parameters called \em omega or \em sigma.
 *
 * \bug This does not affect the start_idx or end_idx! Currently this function
 * is not used, but this would have to be fixed before using it!
 */
void  Parameters::scale(double factor)
{
    for (size_t i = 0; i<size(); i++) {
        if ( (name(i) != "omega") && (name(i) != "sigma") )
            val[i] = val[i].scale(factor);
    }

    //correct();
}

/*!
 * \bug Beware, as this uses idxToD it may underapproximate the result due to
 * rounding errors!
 */
Parameters Parameters::indexExtents() const
{
    Parameters p = *this;
    for (size_t i=0; i<size(); i++) {
        double lower,upper;
        lower = original.idxToD(start_idx[i],i,false);
        upper = original.idxToD(end_idx[i],i,false);
        p.val[i] = Iv(lower,upper);
    }
    return p;
}

double Parameters::volume() const
{
    double vol = 1.0;

    for(size_t i=0; i<val.size(); i++)
        vol *= val[i].width();

    return vol;
}

double Parameters::relativeVolume() const
{
    double vol = 1.0;

    for (size_t i = 0; i<val.size(); i++) {
        vol *= fraction(i);
    }

    return vol;
}

size_t Parameters::id(const std::string &name)
{
    typename std::map<std::string,size_t>::iterator iter = var_id.find(name);
    if (iter == var_id.end()) {
        throw("Parameters::id() - cannot find parameter!");
    }

    return iter->second;
}

std::string Parameters::name(size_t i)
{
    if (i < var_name.size())
        return var_name[i];
    else if (i < (var_name.size() + var_fixed_name.size()) )
        return var_fixed_name[i-var_name.size()];
    else
        throw("Parameters::name - out of bound!");
}

Vec<std::string> Parameters::names()
{
    return var_name;
}

/*!
 * The information includes parameter names, values, initial values, desired
 * resolution, as well as whether the desired resolution has been reached.
 *
 * \todo Implement as stream-operator.
 */
std::string Parameters::print() const
{
    std::string ret;
    ret += "Parameters:\n";
    ret += "Name    | ";
    ret += "Value                      | ";
    ret += "Original Value             | ";
    ret += "Fraction   | ";
    ret += "Minimum Fraction\n";

    for (size_t i = 0; i<size(); i++) {
        char str[100];
        sprintf(str,"%-8s| [ %10.3e, %10.3e ] | [ %10.3e, %10.3e ] | %10.3e | %10.3e",
               var_name[i].c_str(),
               val[i].lower(), val[i].upper(),
               original.val[i].lower(), original.val[i].upper(),
               fraction(i),
               min_fraction[i]);
        ret += str;
        if (minFractionReached(i))
            ret += " * ";
        ret += "\n";
    }
    for (size_t i = 0; i<sizeFixed(); i++) {
        char str[100];
        sprintf(str,"%-8s|   %10.3e                                                                      *\n",
               var_fixed_name[i].c_str(),
               val_fixed[i]);
        ret += str;
    }
    ret += "\n";

    return ret;
}

std::string Parameters::strNames(bool only_sampled) const
{
    std::string out;
    for (size_t i = 0; i<size(); i++) {
        if (only_sampled && Parameters::fraction_set && (min_fraction[i] == 0) )
            continue;
        out.append(var_name[i]).append("; ");
    }
    if (!only_sampled) {
        for (size_t i = 0; i<sizeFixed(); i++) {
            out.append(var_fixed_name[i]).append("; ");
        }
    }
    return out;
}

std::string Parameters::strRanges(bool only_sampled) const
{
    char buffer[256];
    std::string out;
    for (size_t i = 0; i<size(); i++) {
        if (only_sampled && Parameters::fraction_set && (min_fraction[i] == 0) )
            continue;
        sprintf(buffer,"%20.13e, %20.13e; ", val[i].lower(), val[i].upper());
        out.append(buffer);
    }
    if (!only_sampled) {
        for (size_t i = 0; i<sizeFixed(); i++) {
            sprintf(buffer,"%20.13e, %20.13e; ", val_fixed[i], val_fixed[i]);
            out.append(buffer);
        }
    }
    return out;
}

void Parameters::parse(const std::string &str)
{
    ParseParameter parser;
    parser.parse(str);
}

/*!
 * If \a val is an interval of width 0.0 (i.e., a point), then the \a name is
 * added to Parameters::var_fixed_name and \a val to
 * Parameters::original.val_fixed.
 *
 * Otherwise the \a name is added to Parameters::var_name, the \a val to
 * Parameters::original.val, and \a min_fraction to Parameters::min_fraction.
 * Then, the number of quantization steps is determined for the desired relative
 * resolution \a min_fraction. The interval is divided into 2^n steps, where \em
 * n is chosen so that the attained step size is as large as possible while
 * being below the desired requred step size.
 *
 * \bug Parameters::var_id is not updated if a fixed parameter is added, i.e.,
 * if the last parameter that is added is fixed Parameters::var_id is wrong.
 */
void Parameters::add(const std::string &name, const Iv &val, double min_fraction)
{
    if (val.width() == 0.0) {
        Parameters::var_fixed_name.push_back(name);
        Parameters::original.val_fixed.push_back(val.lower());
        return;
    }

    Parameters::var_name.push_back(name);
    Parameters::original.val.push_back(val);
    Parameters::min_fraction.push_back(min_fraction);

    size_t num_steps = 1;
    if (min_fraction > 0.0) {
        Parameters::fraction_set = true;

        while (1.0/(double)num_steps > min_fraction)
            num_steps *= 2;
    }

    Parameters::original.start_idx.push_back(0);
    Parameters::original.end_idx.push_back(num_steps);

    var_id.clear();
    for (size_t i=0; i<var_name.size()+var_fixed_name.size(); i++)
        var_id[Parameters::name(i)] = i;
}

size_t Parameters::dToIdx(double d, size_t dir, bool upper) const
{
    if (!val[dir].includes(d))
        throw("Parameters::dToIdx - value out of bounds!");

    double rel_d = (d-val[dir].lower())/(val[dir].upper()-val[dir].lower());
    double idx_d = (double)start_idx[dir] + rel_d * (double)(end_idx[dir]-start_idx[dir]);

    if (upper)
        return ceil(idx_d);
    else
        return floor(idx_d);
}

/*!
 * \bug This is not necessarily exact and may lead to a smaller than desired
 * parameter-set due to rounding errors.
 */
double Parameters::idxToD(size_t idx, size_t dir, bool upper) const
{
    if ( (idx > end_idx[dir]) || (idx < start_idx[dir]))
        throw("Parameters::idxToD - value out of bounds!");

    if (upper)
        idx++;

    double rel_idx = (double)idx/(double)(end_idx[dir]-start_idx[dir]);
    return (val[dir].lower() + rel_idx*val[dir].width());
}
