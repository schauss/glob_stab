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

#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
#include <map>

#include "iv.h"
#include "vec.hpp"
#include "nd_array.hpp"
#include "taylor_model.h"

/*!
 * \brief The Parameter class represents a set of parameters.
 *
 * This class stores parameter sets, the initial parameter set, the desired
 * resolution per parameter, etc.
 *
 * The initial parameter set as well as the initial parameter values and desired
 * resolution per parameter are stored as static members and set using parse()
 * or add(). Sub-sets of the parameter-set can then be determined using
 * splitAt(), splitAtRelative(), or splitEqual(). Allws checking whether the
 * desired resolution has been reached using minFractionReached(), returns the
 * direction in which the parameter-size in relation with the desired resolution
 * is largest in maxRelFractionDir() and allows to determine the volume() or
 * relativeVolume().
 *
 * \todo Currently one program run can only consider one set of parameters due
 * to the static members used to store parameter names, original parameters,
 * etc. This could also be solved using a shared pointer to the original
 * parameter set or using a separate class (e.g., a sub-class). But this is not
 * the only place we would have to make such changes to get rid of static
 * members.
 */
class Parameters
{
public:
    /*!
     * \brief Create a new parameter set which is set to the initial parameter
     * set given in Parameters::original.
     */
    Parameters();
    /*!
     * \brief Create a new parameter set for which the size is given by \a start
     * and \a end which correspond to steps of the size of the desired
     * resolution.
     * \param start First included step (step size is desired resolution).
     * \param end One after the last included step (step size is desired
     * resolution).
     */
    Parameters(const Vec<size_t> &start, const Vec<size_t> &end);

    //! Return parameter values as Vec of intervals Iv.
    Vec<Iv> value() const;

    //! Return parameter values as as Vec of TaylorModels of type T_coeff.
    template <class T_coeff>
    Vec<TaylorModel<T_coeff> > tmValue() const;


    //! Return parameter values as as Vec of complex TaylorModels of type T_coeff.
    template <class T_coeff>
    Vec<std::complex<TaylorModel<T_coeff> > > tmValueComplex() const;

    //! Return parameters as map from parameter names to intervals Iv.
    std::map<std::string,Iv> map() const;

    //! Return parameters as map from parameter names to TaylorModels of type T_coeff.
    template <class T_coeff>
    std::map<std::string,TaylorModel<T_coeff> > tmMap() const;

    //! Return parameters as map from parameter names to complex TaylorModels of type T_coeff.
    template <class T_coeff>
    std::map<std::string,std::complex<TaylorModel<T_coeff> > > tmMapComplex() const;

    //! Return number of parameters (not including fixed parameters Parameters::val_fixed).
    size_t size() const;
    //! Return number of fixed parameters in Parameters::val_fixed.
    size_t sizeFixed() const {return val_fixed.size();}

    /*!
     * \brief Split parameters along \a direction at \a mid_point.
     * \param direction Direction in which the parameter set should be split,
     * zero-indexed.
     * \param mid_point Absolute value where this parameter should be split.
     * Must be between lower and upper limit of current Parameter value in \a
     * direction.
     * \param middle Better splitting into two equal parts.
     * \return Vec of two Parameters
     */
    Vec<Parameters> splitAt(size_t direction, double mid_point, bool middle=false) const;
    /*!
     * \brief Split parameters along \a direction at relative \a mid_point.
     * \param direction Direction in which the parameter set should be split,
     * zero-indexed.
     * \param mid_point Relative split point between 0.0 (lower bound) and 1.0
     * (upper bound).
     * \return Vec of two Parameters
     */
    Vec<Parameters> splitAtRelative(size_t direction, double mid_point) const;
    /*!
     * \brief Split parameters along \a direction into \a num_parts parts.
     * \param direction Direction in which the parameter set should be split,
     * zero-indexed.
     * \param num_parts How many parts should the parameter set be split into?
     * \return Vec of \a num_parts Parameters
     */
    Vec<Parameters> splitEqual(size_t direction, size_t num_parts=2) const;

    //! Return intersection with Parameters \a p as Parameters.
    Parameters intersect(const Parameters &p) const;

    //! Return width of interval in direction \a i (relative to initial width).
    double  fraction(size_t i) const;
    //! Return desired resolution in direction \a i (relative to initial width).
    double  minFraction(size_t i) const;
    /*!
     * \brief Has the desired resolution been reached in direction \a i (or in
     * all directions if \a i<0, i.e., if \a i not given.)
     */
    bool    minFractionReached(int i = -1) const;

    //! Returns direction of largest parameter (relative to desired resolution).
    size_t  maxRelFractionDir() const;

    //! True if parameter in direction \a i is identical to initial interval.
    bool    isOriginal(size_t i) const;

    //! Intersect parameter set with initial parameter set.
    void    correct();
    /*!
     * \brief Scale parameter set by a \a factor, i.e., mid-points of each
     * parameter stay the same but the intervals are scaled around this mid
     * point.
     */
    void    scale(double factor);

    /*!
     * \brief Return parameter set corresponding to current index.
     *
     * This should always enlarges the parameter set. It is used to determine
     * for which region we should store zero-inclusion in OptimizationResult.
     */
    Parameters indexExtents() const;

    //! Volume of parameter set (product of interval widths).
    double volume() const;
    //! Relative volume of parameter set (product of fractions, see fraction()).
    double relativeVolume() const;

    //! Return information on parameter set as string.
    std::string print() const;
    /*!
     * \brief  Return all parameter names as string (separated by semi-colons).
     * \param only_sampled If true, then only interval parameters with a desired
     * resolution are returned. Otherwise all interval parameters and fixed
     * parameters are returned.
     */
    std::string strNames(bool only_sampled=false) const;
    /*!
     * \brief  Return all parameter ranges as string (separated by semi-colons).
     * \param only_sampled If true, then only interval-parameters with a desired
     * resolution are returned. Otherwise all interval parameters and fixed
     * parameters are returned.
     */
    std::string strRanges(bool only_sampled=false) const;

    /*!
     * \brief Parse string \a str by calling ParseParameter::parse().
     *
     * After parsing the string the function ParseParameter::parse() adds the
     * parameter(s) to the initial parameter set by calling Parameters::add().
     * See the ParseParameter class for more details on required syntax, etc.
     */
    static void parse(const std::string &str);
    /*!
     * \brief Add a parameter to the initial parameter set.
     * \param name Parameter name
     * \param val Parameter value
     * \param min_fraction Desired resolution for boundary mapping (relative to
     * width of initial interval given in \a val.
     */
    static void add(const std::string &name, const Iv &val, double min_fraction = 0.0);

    //! Index of parameter with \a name (works for intervals and fixed).
    static size_t id(const std::string &name);
    //! Name of parameter with index \a i (works for intervals and fixed).
    static std::string name(size_t i);
    //! Vec of all parameter names (intervals and fixed).
    static Vec<std::string> names();

    // Start index of parameter set, same dimension as number of interval parameters.
    const Vec<size_t> &startIdx() const {return start_idx;}
    // End index of parameter set, same dimension as number of interval parameters.
    const Vec<size_t> &endIdx() const {return end_idx;}

private:
    /*!
     * \brief Private constructor to be called to create static initial
     * parameters Parameters::original.
     * \param stat Value not relevant, only important to call with bool-type.
     */
    Parameters(bool stat);

    /*!
     * \brief Return index of quantized paramter for given value.
     * \param d Parameter value for which index should be found.
     * \param dir Of which parameter should the index be determined.
     * \param upper Increment returned value by one.
     */
    size_t dToIdx(double d, size_t dir, bool upper = false) const;
    /*!
     * \brief Return parameter value of quantized parameter for given index.
     * \param idx Index for which the paramter value should be returned.
     * \param dir Of which parameter should the value be determined.
     * \param upper Value of upper bound of quantization step \a idx? Otherwise
     * return lower bound!
     */
    double idxToD(size_t idx, size_t dir, bool upper = false) const;

    //! Interval values of parameter set. Names given in Parameters::var_names.
    Vec<Iv> val;
    //! Fixed values (double) of parameter set. Names given in Parameters::var_fixed_names.
    Vec<double> val_fixed;
    //! Start index of quantized parameters for this parameter set.
    Vec<size_t> start_idx;
    //! End index of quantized parameters for this parameter set (one after end)
    Vec<size_t> end_idx;

    //! Initial parameter set. Initialized using parse() or add()
    static Parameters original;
    //! Names of interval parameters in Parameters::val.
    static Vec<std::string> var_name;
    //! Names of fixed parameters in Parameters::val_fixed.
    static Vec<std::string> var_fixed_name;
    //! Map from parameter names (interval+fixed) to parameter index.
    static std::map<std::string,size_t> var_id;

    //! Desired resolution (relative to initial interval) for Parameters::val.
    static Vec<double> min_fraction;
    //! Has a desired resolution been specified for any variable?
    static bool fraction_set;
};

#endif // PARAMETER_H
