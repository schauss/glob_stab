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

#include "bernstein.h"
#include "polynomial.h"
#include "taylor_model.h"
#include "helpers.h"

#include <cstdlib>

template <class T>
size_t Bernstein<T>::count = 0;

template <class T>
size_t Bernstein<T>::max_count = 0;

template <class T>
Bernstein<T>::Bernstein() :
    coeff_done(false),
    iv_bound_done(false),
    empty(true),
    ex(BernsteinExponent::getInstance(std::vector<uint8_t>(1,0)))
{
}

/*!
 * BEWARE: This constructor assumes the Polynomial is on the domain [0,1] and does not actually ensure this!
 */
template <class T>
Bernstein<T>::Bernstein(const Polynomial<T> &pol, std::vector<uint8_t> degree) :
    coeff_done(false),
    iv_bound_done(false),
    empty(false),
    ex(BernsteinExponent::getInstance(degree.empty()?pol.maxOrders():degree))
{
    poly = Vec<T>(ex->maxId(),0.0);
    for (typename std::set<Monomial<T> >::const_iterator it = pol.val().begin(); it!=pol.val().end(); it++) {
        size_t idx = ex->toId(it->exp());
        poly[idx] = it->coeff() * ex->binom_vec[idx];
    }

    count++;
    if (count > max_count)
        max_count = count;
}

/*!
 * This constructor is needed for the subdivision-algorithm implemented in splitAt().
 */
template <class T>
Bernstein<T>::Bernstein(const Vec<T> pol, BernsteinExponent *exp) :
    coeff_done(true),
    iv_bound_done(false),
    empty(false),
    ex(exp)
{
    poly = pol;

    count++;
    if (count > max_count)
        max_count = count;
}

template <class T>
Bernstein<T>::Bernstein(const Bernstein<T> &b) :
    poly(b.poly),
    coeff_done(b.coeff_done),
    iv_bound(b.iv_bound),
    iv_inner_bound(b.iv_inner_bound),
    iv_bound_done(b.iv_bound_done),
    empty(b.empty),
    derivative(b.derivative),
    max_derivative_dir(b.max_derivative_dir),
    ex(b.ex)
{
    if (!empty) {
        count++;
        if (count > max_count)
            max_count = count;
    }
}

template <class T>
Bernstein<T>::~Bernstein()
{
    if (!empty)
        count--;
}

template <class T>
Bernstein<T>& Bernstein<T>::operator=(const Bernstein<T> &b)
{
    poly                = b.poly;
    coeff_done          = b.coeff_done;
    iv_bound            = b.iv_bound;
    iv_inner_bound      = b.iv_inner_bound;
    iv_bound_done       = b.iv_bound_done;
    derivative          = b.derivative;
    max_derivative_dir  = b.max_derivative_dir;
    ex                  = b.ex;

    if (empty && !b.empty) {
        count++;
        if (count > max_count)
            max_count = count;
    }
    else if (!empty && b.empty)
        count--;

    empty = b.empty;

    return *this;
}

/*!
 * This is determined by comparing the outer approximation and inner approximation.
 * If the upper bounds are identical, then the upper bound returned by bound() or innerBound() is exact.
 */
template <class T>
bool Bernstein<T>::upperSharp() const
{
    if (isEmpty())
        throw("Bernstein::upperSharp() called on empty Bernstein polynomial!");

    if (!iv_bound_done)
        calcBound();

    return ( bound().upper() == innerBound().upper() );
}

/*!
 * This is determined by comparing the outer approximation and inner approximation.
 * If the lower bounds are identical, then the lower bound returned by bound() or innerBound() is exact.
 */
template <class T>
bool Bernstein<T>::lowerSharp() const
{
    if (isEmpty())
        throw("Bernstein::lowerSharp() called on empty Bernstein polynomial!");

    if (!iv_bound_done)
        calcBound();

    return ( bound().lower() == innerBound().lower() );
}

/*!
 * The outer approximation is calculated using calcBound() if this has not been done yet.
 */
template <class T>
Iv Bernstein<T>::bound() const
{
    if (isEmpty())
        throw("Bernstein::bound() called on empty Bernstein polynomial!");

    if (!iv_bound_done)
        calcBound();

    return iv_bound;
}

/*!
 * The inner approximation is calculated using calcBound() if this has not been done yet.
 */
template <class T>
Iv Bernstein<T>::innerBound() const
{
    if (isEmpty())
        throw("Bernstein::innerBound() called on empty Bernstein polynomial!");

    if (!iv_bound_done)
        calcBound();

    return iv_inner_bound;
}

/*!
 * The range for each dimension is calculated using calcRange if this has not been done yet.
 */
template <class T>
Vec<double> Bernstein<T>::range() const
{
    if (isEmpty())
        throw("Bernstein::range() called on empty Bernstein polynomial!");

    if (directional_range.empty())
        calcRange();

    return directional_range;
}

/*!
 * The partial derivative is calculated using calcDerivative() if this has not been done yet.
 */
template <class T>
Vec<double> Bernstein<T>::dx() const
{
    if (isEmpty())
        throw("Bernstein::dx() called on empty Bernstein polynomial!");

    if (derivative.empty())
        calcDerivative();

    return derivative;
}

/*!
 * These Bernstein coefficients are exact, i.e., for the corresponding variable values (maximum or minimum for each
 * dimension) the polynomial has the value of the Bernstein coefficient.
 */
template <class T>
Vec<T> Bernstein<T>::coeffs_vertex() const
{
    bound();
    Vec<T> ret;
    for (size_t i=0; i<poly.size(); i++) {
        if (ex->isVertex(i))
            ret.push_back(poly[i]);
    }
    return ret;
}

/*!
 * The partial derivative is calculated using calcDerivative() if this has not been done yet.
 */
template <class T>
int Bernstein<T>::maxDxDir() const
{
    if (isEmpty())
        throw("Bernstein::maxDxDir() called on empty Bernstein polynomial!");

    if (derivative.empty())
        calcDerivative();

    return max_derivative_dir;
}

template <class T>
Vec<uint8_t> Bernstein<T>::getOrders() const
{
    return ex->max_degree;
}

template <class T>
unsigned int Bernstein<T>::getIdx(Vec<uint8_t> &exp) const
{
    if (exp.size() != ex->numVar())
        throw("Bernstein::getIdx called with wrong length of exponent-vector!");

    uint8_t exponents[ex->numVar()];
    for (size_t i=0; i<ex->numVar(); i++)
        exponents[i] = exp[i];

    return (ex->toId(exponents));
}

/*!
 * This algorithm calculates the Bernstein coefficients from the coefficients of a polynomial on the unit interval box
 * efficiently. It is taken from \cite Garloff1986.
 */
template <class T>
void Bernstein<T>::fastCoeff() const
{
    if (coeff_done)
        throw("Bernstein::fastCoeff - Already done!");

    for (int dim = 0; dim < (int)ex->numVar(); dim++) {
        int i = 0;
        int step = ex->step[dim];
        int order = ex->max_degree[dim];
        while (i < (int)ex->maxId())
        {
            for (int j=0; j<step; j++)  {
                for (int k_start=0; k_start < order; k_start++) {
                    for (int k=order-1; k>=k_start; k--) {
                        int ind = i+j+k*step;
                        poly[ind+step] += poly[ind];
                    }
                }
            }
            i += step*(order+1);
        }
    }

    coeff_done = true;
}

/*!
 * This algorithm calculates the Bernstein coefficients of two subboxes when splitting along dimension dim at the point
 * split_point.
 * It is taken from \cite Zettler1998.
 */
template <class T>
Vec<Bernstein<T> > Bernstein<T>::splitAt(size_t dim, double split_point) const
{
    if (isEmpty())
        throw("Bernstein::subdivide() called on empty Bernstein polynomial!");

    if (!coeff_done)
        throw("Bernstein::subdivide - Coeffs not done!");

    Vec<Bernstein<T> > ret;

    Vec<T> poly_A = poly;
    Vec<T> poly_B(ex->maxId(),0.0);

    int i = 0;
    int step = ex->step[dim];
    int order = ex->max_degree[dim];
    while (i < (int)ex->maxId())
    {
        for (int j=0; j<step; j++)  {
            // extract first coefficient for poly_B here (order i_r==n_r)
            poly_B[i+j+order*step] = poly_A[i+j+order*step];

            // k_start+1 corresponds to k in Garloff-Paper
            for (int k_start=0; k_start < order; k_start++) {
                // inner loop makes additional storage matrices unnecessary (compute from back to front)
                for (int k=order-1; k>=k_start; k--) {
                    int ind = i+j+k*step;
                    poly_A[ind+step] = (1-split_point)*poly_A[ind]+ split_point*poly_A[ind+step];
                }
                // here all values for one k (Garloff-Paper are computed and we can extract the values for poly_B)
                // k_start-1 as his k started at 1!
                int ind = i+j+(order-k_start-1)*step;
                poly_B[ind] = poly_A[i+j+order*step];
            }
        }
        i += step*(order+1);
    }

    ret.push_back(Bernstein<T>(poly_A,ex));
    ret.push_back(Bernstein<T>(poly_B,ex));
    return ret;
}

/*!
 * This algorithm estimates the maximum partial derivative allong each dimension.
 * It is taken from \cite Zettler1998.
 */
template <class T>
void Bernstein<T>::calcDerivative() const
{
    if(!derivative.empty())
        throw("Bernstein::calcDerivatives - already calculated!");

    derivative = Vec<double>(ex->numVar(),0.0);
    double overall_max_derivative = 0.0;

    for (int dim = 0; dim < (int)ex->numVar(); dim++) {
        double max_derivative = 0.0;
        int i = 0;
        int step = ex->step[dim];
        int order = ex->max_degree[dim];
        while (i < (int)ex->maxId())
        {
            for (int j=0; j<step; j++)  {
                for (int k=0; k < order; k++) {
                    int ind = i+j+k*step;
                    double delta = deltaMax(ind+step,ind);
                    if (delta > max_derivative) {
                        max_derivative = delta;
                    }
                }
            }
            i += step*(order+1);
        }
        derivative[dim] = order * max_derivative;
        if (derivative[dim] > overall_max_derivative) {
            overall_max_derivative = derivative[dim];
            max_derivative_dir = dim;
        }
    }
}

/*!
 * This algorithm calculates the maximum difference between minimum and maximum Bernstein coefficient for each dimension
 * while the indices of all other dimensions are kept constant.
 */
template <class T>
void Bernstein<T>::calcRange() const
{
    if(!directional_range.empty())
        throw("Bernstein::calcRange - already calculated!");

    directional_range = Vec<double>(ex->numVar(),0.0);

    for (int dim = 0; dim < (int)ex->numVar(); dim++) {
        double max_range =  -INFINITY;

        size_t i = 0;
        size_t step = ex->step[dim];
        size_t order = ex->max_degree[dim];
        while (i < ex->maxId())
        {
            for (size_t j=0; j<step; j++)  {
                double min_val = INFINITY;
                double max_val = -INFINITY;
                for (size_t k=0; k <= order; k++) {
                    size_t ind = i+j+k*step;
                    if (min_val > minVal(ind))
                        min_val = minVal(ind);
                    if (max_val < maxVal(ind))
                        max_val = maxVal(ind);
                }
                if (max_range < (max_val - min_val))
                    max_range = (max_val - min_val);
            }
            i += step*(order+1);
        }
        directional_range[dim] = max_range;
    }
}

/*!
 * The outer approximation of the range of the polynomial is given by the minimum and maximum Bernstein coefficient.
 * An inner approximation of the range of the polynomial is given by the minimum and maximum Bernstein coefficient of
 * the vertices of the interval-box (i.e., the minimum and maximum value for each parameter).
 * For details see e.g., \cite Zettler1998.
 */
template <class T>
void Bernstein<T>::calcBound() const
{
    if (iv_bound_done)
        throw("Bernstein::calcBounds - Bounds already calculated!");

    if (!coeff_done)
        fastCoeff();

    double lower = minVal(0);
    double upper = maxVal(0);
    for (size_t i = 1; i<poly.size(); i++) {
        if (minVal(i) < lower) {
            lower = minVal(i);
        }
        if (maxVal(i) > upper) {
            upper = maxVal(i);
        }
    }
    iv_bound = Iv(lower,upper);

    Vec<size_t> vertex_id = ex->vertexId();
    double inner_lower = maxVal(vertex_id[0]);
    double inner_upper = minVal(vertex_id[0]);
    for (size_t i = 1; i<vertex_id.size(); i++) {
        size_t idx = vertex_id[i];
        if (maxVal(idx) < inner_lower) {
            inner_lower = maxVal(idx);
        }
        if (minVal(idx) > inner_upper) {
            inner_upper = minVal(idx);
        }
    }
    iv_inner_bound = Iv(inner_lower,inner_upper);

    iv_bound_done = true;
}

/*!
 * Generally this simply returns the absolute value of the difference between the values of the two Bernstein
 * coefficients.
 */
template <class T>
double Bernstein<T>::deltaMax(size_t id1, size_t id2) const
{
    return fabs(poly[id1]-poly[id2]);
}

/*!
 * In case the the Bernstein coefficients are intervals this function returns the maximum absolute difference between
 * any value of the two intervals.
 */
template <>
double Bernstein<Iv>::deltaMax(size_t id1, size_t id2) const
{
    Iv diff = poly[id1]-poly[id2];
    double tmp1 = fabs(diff.lower());
    double tmp2 = fabs(diff.upper());

    return (tmp1>tmp2?tmp1:tmp2);
}

/*!
 * Generally this simply returns the value of the Bernstein coefficient.
 */
template <class T>
double Bernstein<T>::maxVal(size_t id) const
{
    return poly[id];
}

/*!
 * In case the Bernstein coefficients are intervals this function returns the upper bound of the coefficient interval.
 */
template <>
double Bernstein<Iv>::maxVal(size_t id) const
{
    return poly[id].upper();
}

/*!
 * Generally this simply returns the value of the Bernstein coefficient.
 */
template <class T>
double Bernstein<T>::minVal(size_t id) const
{
    return poly[id];
}

/*!
 * In case the Bernstein coefficients are intervals this function returns the lower bound of the coefficient interval.
 */
template <>
double Bernstein<Iv>::minVal(size_t id) const
{
    return poly[id].lower();
}

template class Bernstein<double>;
template class Bernstein<Iv>;

template <class T>
std::ostream & operator<<(std::ostream &os, const Bernstein<T>& b)
{
    os << "bound: " << b.bound() << std::endl;
    os << "dx:    " << b.dx() << std::endl;
    os << "range: " << b.range() << std::endl;
    os << "order: " << Vec<int>(b.getOrders()) << std::endl;
    os << "coefficients: " << std::endl;
    for (size_t i=0; i<b.poly.size(); i++) {
        std::cout << Vec<int>(b.ex->toVec(i)) << "\t" << b.poly[i] << std::endl;
    }
    return os;
}

template std::ostream & operator<<(std::ostream &os, const Bernstein<double>& b);
template std::ostream & operator<<(std::ostream &os, const Bernstein<Iv>& b);
