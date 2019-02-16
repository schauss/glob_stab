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

#ifndef FUNCTION_H
#define FUNCTION_H

#include <boost/shared_ptr.hpp>

#include "bernstein.h"
#include "bernstein_2d.h"
#include "equation.h"
#include "parameters.h"
#include "taylor_model.h"
#include "value_set.h"

/*!
 * \brief The class Function can either be used to evaluate an Equation or a ValueSet.
 *
 * It stores the Equation / ValueSet as well as a set of parameters and offers a common interface to zero-inclusion and
 * -exclusion check, maximum partial derivative, partial derivative of the size of the interval remainder.
 *
 */
class Function
{
public:
    /*!
     * \brief Construct Function from Equation and Parameters
     * \param eq Equation to be wrapped by function
     * \param param Parameters for which the Equation is evaluated
     *
     * Does not actually evaluate the Equation!
     */
    Function(const Equation &eq, const Parameters &param = Parameters());
    /*!
     * \brief Construct Function from ValueSet and Parameters
     * \param vs ValueSet to be wrapped by function
     * \param param Parameters for which the ValueSet is evaluated
     *
     * Does not actually evaluate the ValueSet!
     */
    Function(const ValueSet &vs, const Parameters &param = Parameters());
    //! Destructor.
    virtual ~Function() {}

    /*!
     * \brief Evaluate Taylor Model of Equation or ValueSet
     *
     * This results in Bernstein polynomials and an interval remainder.
     *
     * Does nothing if Function has already been evaluated_from_euqation!
     */
    void evaluateFromEquation();
    /*!
     * \brief Evaluate the Equation or ValueSet.
     *
     * This means that inner and outer approximations are determined and zero-inclusion and -exclusion are evaluated.
     */
    void evaluate();
    //! Clear Bernstein coefficients in poly or poly2d and set function to be not-evaluated.
    void clear();

    //! Return vector of values for Equation, not implemented for ValueSet
    Vec<Iv> value() const;
    //! Return parameters
    const Parameters& parameters() const {return param;}

    /*!
     * \brief Return estimated maximum partial derivative of the polynomial part of the ValueSet or Equation
     *
     * Maximum partial derivative means that over the complete domain of the Bernstein polynomial this is the estimated
     * maximum of the partial derivative for each dimension. In case of a ValueSet the square norm of the two estimated
     * maximum partial derivatives is returned (i.e., the square of the absolute value of the complex value)
     */
    Vec<Vec<double> > dpoly_dparam() const;
    /*!
     * \brief Return estimated maximum partial derivative of the interval remainder of the ValueSet or Equation
     *
     * The parameter set is split once in each direction and the Taylor Model of the Equation / Value Set is evaluated.
     * Then the size of the interval remainder is returned. This is computationally very expensive!
     *
     * Note that in case of the ValueSet the square norm of the size of the complex interval remainder is returned.
     *
     * \todo TIn case of an Equation this function returns the larger of the two interval remainder widths. In case of a
     * ValueSet we actually only evaluate one Taylor Model in each dimension (the lower half) to reduce computation
     * time. Obviously, this is not ideal!
     */
    Vec<Vec<double> > drest_dparam() const;

    //! Has the function been evaluated (has evaluate() been called)?
    bool isEvaluated() const {return evaluated;}
    //! Has the function been evaluated from the Equation or ValueSet (has evaluateFromEquation() been called)?
    bool isEvaluatedFromEquation() const {return evaluated_from_equation;}

    /*!
     * \brief Check whether the function include zero, i.e., is the negative interval remainder included in an inner
     * approximation.
     *
     * In case of an Equation this indicates that "one" of the Equations (it is a vector) has this property.
     */
    bool includesZero() const;
    /*!
     * \brief Check whether the negative interval remainder is included in the outer approximation of the polynomial
     * part but not in an inner approximation.
     *
     * In case of an Equation this indicates that "one" of the Equations (it is a vector) has this property!
     */
    bool partiallyIncludesZero1() const;
    /*!
     * \brief Check whether zero is included in the outer approximation of the Function.
     *
     * In case of an Equation this indicates that "one" of the Equations (it is a vector) has this property!
     */
    bool partiallyIncludesZero2() const;
    /*!
     * \brief Checks whether zero is excluded from the outer approximation of the Function, i.e., the negative interval
     * remainder is excluded from the outer approximation of the polynomial part.
     *
     * \todo In case of an Equation this looks strange: it checks whether any of the Equations is less or equal to zero
     * or all of the Equations are larger than zero ???
     */
    bool excludesZero() const;

    /*!
     * \brief Checks whether the boundary of an inner approximation or outer approximation of the polynomial part or
     * of the outer approximation of the overall Function is zero.
     *
     * In case of an Equation this indicates that "one" of the Equations (it is a vector) has this property!
     *
     * \todo I have no idea where the use-case of this function is.
     * \bug Currently this is not implemented for ValueSets and always returns false!
     */
    bool boundaryIsZero() const;
    //! Returns true if all remainder intervals are finit.
    bool restIsFinite() const;

    //! Are we evaluating a ValueSet (true) or an Equation (false)
    bool doValueSet() const {return (value_set.get()!=0);}
    //! Are we evaluating the boundary of a ValueSet (true) or not (false)
    bool doValueSetBoundary() const {return ( (value_set.get()!=0) && (value_set->doBoundary()) );}

    /*!
     * \brief Split Function at given mid_point.
     * \param direction Dimrection in which the Function (or better the Parameter set) is split.
     * \param mid_point Relative point along the dimension where the Parameter set is split (defaults to the midpoint).
     *
     * This returns two Functions. The interval remainder is set to the interval remainder of this function while a
     * Bernstein Subdivision Algorithm is used to determine the Bernstein coefficients after the split.
     */
    std::vector<boost::shared_ptr<Function> > splitAt(size_t direction, double mid_point = 0.5) const;
    //! Clear the Bernstein coefficients and slightly enlarge the parameter set
    void enlarge();

    //! Check value of has_marginal
    bool hasMarginal();
    //! Set value of has_marginal
    void setMarginal(bool has_marginal);

private:
    //! Vector of Bernstein polynomials with interval coefficients. Either this or poly2d exist.
    Vec<Bernstein<Iv> >   poly;
    //! Vector of Bernstein2d polynomials (complex Bernstein polynomials) with interval coefficients. Either this or poly exist.
    Vec<Bernstein2d<Iv> > poly2d;
    //! Interval remainder
    Vec<Iv>               rest;
    //! Vector specifying for which Polynomial in poly 0 is excluded from the outer approximation.
    Vec<bool>             done;

    //! Equation to evaluate. Results in poly. Either this or value_set is given.
    boost::shared_ptr<const Equation>   equation;
    //! Value set to evaluate. Results in poly2d. Either this or equation is given.
    boost::shared_ptr<const ValueSet>   value_set;
    //! Parameters for which the equation or value set should be evaluated.
    Parameters                    param;

    //! Function has been evaluated (i.e., bounds have been determined)
    bool evaluated;
    //! Function has been evaluated from equation or value set (i.e., interval remainder was determined for this set of parameters, not for bigger set)
    bool evaluated_from_equation;
    //! Function is marginally stable. \todo Not nice to store this here!
    bool has_marginal;

    //! Stream operator for debug-output
    friend std::ostream & operator<<(std::ostream &os, Function& fcn);
};

//! Stream operator for debug-output
std::ostream & operator<<(std::ostream &os, Function& fcn);

#endif // FUNCTION_H
