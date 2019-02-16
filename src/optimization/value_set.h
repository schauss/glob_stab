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

#ifndef VALUE_SET_H
#define VALUE_SET_H

#include <complex>
#include <map>
#include <ostream>
#include <string>

#include "equation.h"
#include "vec.hpp"

/*!
 * \brief The class ValueSet evaluates the complex value set of an equation.
 *
 * The equation The class makes use of ParseEquation to evaluate the equation on calls to value() or valueMat().
 */
class ValueSet
{
public:
    /*!
     * \brief Construct ValueSet and initialize member variables.
     * \param chEq Set class member ch_eq to chEq. This is a vector of coefficients and is evaluated as
     * ch_eq[0]+ch_eq[1]*s+ch_eq[2]*s^2+...
     * \param delta Desired damping. Note that damping!=0 is problematic for systems with time delay if the inverse
     * value set is evaluated.
     *
     * Does not perform any actual evaluation.
     */
    ValueSet(const Equation &chEq, double delta);
    //! Destructor.
    virtual ~ValueSet() {}

    /*!
     * \brief Evaluate value set for given parameters and return result.
     * \param param Map of uncertain parameters
     * \return A Vec<T> with two elements where the first is the real part and the second the imaginary part of the
     * value set.
     */
    template <class T>
    Vec<T> value(const std::map<std::string,T> &param) const;

    /*!
     * \brief Evaluate value set of ch_eq[0] for s=0.
     * \param param Map of uncertain parameters
     * \return A real value T (as s=0 the result must be real)
     */
    template <class T>
    T valueS0(const std::map<std::string,T> &param) const;

    bool doBoundary() const {return boundary;}  //!< Evaluate stability crossing set, i.e., set sigma=0 ?
    void doBoundary(bool val) {boundary = val;} //!< Evaluate stability crossing set, i.e., set sigma=0 if val==True !

    bool doComplex() const {return complex;}    //!< Set s=sigma+j*omega (necessary for systems with time delay) ?
    void doComplex(bool val) {complex = val;}   //!< Set s=sigma+j*omega if val==True !

    bool doInverse() const {return inverse;}    //!< Evaluate inverse value set (allows evaluating unbounded value set) ?
    void doInverse(bool val) {inverse = val;}   //!< Evaluate inverse value set if val==True !

private:
    /*!
     * \brief Retrieve parameter value from parameter map
     * \param param Parameter map
     * \param name Parameter name
     * \return Parametere value
     */
    template <class T>
    T get(const std::map<std::string,T> &param, const std::string &name) const;

    /*!
     * \brief Identifiers "s" in the coefficients of the characteristic equation are replaced by this.
     * \param sigma Real part of s (not yet modified according to delta-value, is done within function)
     * \param omega Imaginary part of s
     * \return Complex value s (or s^(-1)) corrected by delta-shift
     *
     * This replacement is only performed if doComplex()==true. Depending on inverse and the value of delta different
     * results are returned. Note, that a value of delta!=0.0 together with an inverse results in an infinite s!
     */
    template <class T>
    std::complex<T> s(const T &sigma, const T &omega) const;

    /*!
     * \brief Returns powers of s. The coefficients of the characteristic equation are multiplied by these.
     * \param sigma Real part of s (actually sigma-delta*omega is expected, i.e., delta is not taken into account within
     * this function).
     * \param omega Imaginary part of s
     * \param i Index of coefficient
     * \param n Maximum power of s, i.e., number of coefficients in ch_eq - 1
     * \return Complex value s^i (doInverse()==false) or s^(n-i) (doInverse()==true)
     *
     * The expression (sigma + j*omega)^i or (signa + j*omega)^(n-i) is evaluated. We explicitly implement the binomial
     * expression as this results in tighter bounds than when taking the complex s-value to some power!
     */
    template <class T>
    std::complex<T> s(const T &sigma, const T &omega, size_t i, size_t n) const;

    Equation ch_eq; //!< Coefficients of characteristic equation

    bool boundary;  //!< Evaluate boundary (i.e., set sigma=0 or sigma = -omega*delta)
    bool complex;   //!< Set s=sigma+j*omega when evaluating the coefficients of the characteristic equation
    bool inverse;   //!< Evaluate the inverse value set

    double delta;   //!< Value of delta for desired damping

    //! Stream operator for debug-output
    friend std::ostream & operator<<(std::ostream &os, const ValueSet& vs);
};

//! Stream operator for debug-output
std::ostream & operator<<(std::ostream &os, const ValueSet& vs);

#endif // VALUE_SET_H
