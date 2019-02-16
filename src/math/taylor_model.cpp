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

#include "taylor_model.h"
#include "bernstein.h"

template <class T>
Iv TaylorModel<T>::domain(0,1);
template <class T>
int TaylorModel<T>::max_total_degree = 10;
template <class T>
bool TaylorModel<T>::correct_trig = false;

template <class T>
boost::program_options::options_description TaylorModel<T>::options()
{
    // Declare generic options.
    boost::program_options::options_description opt("Taylor-Model Options");
    opt.add_options()
        ("tm_order",     boost::program_options::value<int>()->default_value(10), "order of taylor model")
        ("tm_center", "use domain [-1,1] instead of [0,1]")
        ("tm_b_off", "turn off bernstein bounding")
    ;
    return opt;
}

template <class T>
void TaylorModel<T>::init(const Vec<std::string> &var_names)
{
    max_total_degree = CONFIG_VAR(tm_order,int);

    if (CONFIG_ISSET(tm_center))
        domain = Iv(-1,1);
    else
        domain = Iv(0,1);

    Polynomial<T>::doBernsteinBound(!CONFIG_ISSET(tm_b_off));

    std::cout << "Using TaylorModel of degree " << max_total_degree << " on domain " << domain << (Polynomial<double>::doBernsteinBound()?" with Bernstein bounding!":" with Interval bounding!") << std::endl;

    std::cout << "Setting up " << var_names.size() << " variables for polynomial!" << std::endl;

    Monomial<T>::fixVariables(var_names);

    std::cout << std::endl;
}

template <class T>
TaylorModel<T>::TaylorModel(bool not_exactly_zero)
    : rest(false)
{
    if (not_exactly_zero)
        throw("TaylorModel::TaylorModel(true) not valid!");

}

template <class T>
TaylorModel<T>::TaylorModel(const std::string& var)
    : poly(var), rest(false)
{
}

template <class T>
TaylorModel<T>::TaylorModel(const char var[])
    : poly(var), rest(false)
{
}

template <class T>
TaylorModel<T>::TaylorModel(double val)
    : poly(val), rest(false)
{
}

template <class T>
TaylorModel<T>::TaylorModel(const Iv &iv)
    : rest(iv)
{
}

/*!
 * \bug May underapproximate due to possible rounding errors.
 */
template <class T>
TaylorModel<T>::TaylorModel(const std::string& var, const Iv &i)
    : rest(false)
{
    double offset;
    double scale;

    if (domain == Iv(-1,1)) {
        offset = 0.5*(i.upper()+i.lower());
        scale = 0.5*(i.upper()-i.lower());
    }
    else if (domain == Iv(0,1)){
        offset = i.lower();
        scale = i.upper()-i.lower();
    } else
        throw("TaylorModel::TaylorModel(std::string,Iv) - unsupported domain!");

    poly = Polynomial<T>(offset) + Polynomial<T>(scale) * Polynomial<T>(var);
}

template <class T>
TaylorModel<T>::~TaylorModel()
{

}

template <class T>
const Polynomial<T> &TaylorModel<T>::polynomial_01() const
{
    if (domain == Iv(0,1))
        return poly;

    if (poly_01.isEmpty())
        poly_01 = poly.changeDomain(domain,Iv(0,1));

    return poly_01;
}

template <class T>
Iv TaylorModel<T>::bound() const
{
    return poly.bound(domain)+rest;
}

/*!
 * This operation adds the rhs-Polynomial to this Polynomial and the
 * rhs-interval-remainder to this interval remainder.
 */
template <class T>
TaylorModel<T> &TaylorModel<T>::operator+= (const TaylorModel<T> &rhs)
{
    poly += rhs.poly;
    rest += rhs.rest;

    if (!poly_01.isEmpty())
        poly_01 = Polynomial<T>();

    return *this;
}

/*!
 * Negates the rhs-TaylorModel using operator-() and then adds it to this
 * TaylorModel using operator+=().
 */
template <class T>
TaylorModel<T> &TaylorModel<T>::operator-= (const TaylorModel<T> &rhs)
{
    return (*this += -rhs);
}

/*!
 * Multiplies the two Polynomials and limits the maximum degree of the exponents
 * to the value given in max_total_degree using cutoff(). The Polynomial part
 * of higher degree is part of the inerval remainder. In addition, the interval
 * remainder consists of the product of the two interval remainders with each
 * as well as with an outer approximation of the Polynomial of the respective
 * other TaylorModel.
 */
template <class T>
TaylorModel<T> &TaylorModel<T>::operator*= (const TaylorModel<T> &rhs)
{
    rest = (rest*rhs.poly.bound(domain))+(rhs.rest*poly.bound(domain)) + rest*rhs.rest;

    poly *= rhs.poly;
    cutoff();

    if (!poly_01.isEmpty())
        poly_01 = Polynomial<T>();

    return *this;
}

/*!
 * Multiplies the TaylorModel y the reciprocal of the rhs-TaylorModel which is
 * determined using inv().
 */
template <class T>
TaylorModel<T> &TaylorModel<T>::operator/= (const TaylorModel<T> &rhs)
{
    return (*this *= rhs.inv());
}

template <class T>
TaylorModel<T> &TaylorModel<T>::operator+= (double rhs)
{
    *this += TaylorModel<T>(rhs);
    return *this;
}

template <class T>
TaylorModel<T> &TaylorModel<T>::operator-= (double rhs)
{
    *this -= TaylorModel<T>(rhs);
    return *this;
}

template <class T>
TaylorModel<T> &TaylorModel<T>::operator*= (double rhs)
{
    *this *= TaylorModel<T>(rhs);
    return *this;
}

template <class T>
TaylorModel<T> &TaylorModel<T>::operator/= (double rhs)
{
    *this /= TaylorModel<T>(rhs);
    return *this;
}

template <class T>
TaylorModel<T> TaylorModel<T>::operator-() const
{
    TaylorModel ret;
    ret.poly = -poly;
    ret.rest = -rest;

    return ret;
}

/*!
 * The reciprocal can be determined using the horner-scheme implemented in
 * horner(). For more details, see \cite Makino2003.
 */
template <class T>
TaylorModel<T> TaylorModel<T>::inv() const
{
    if (max_total_degree<0)
        throw("TaylorModel::Inv - cannot be used with taylor-model of unlimited total degree!");

    Iv c_f = cF().mid();
    TaylorModel<T> f_bar = *this - c_f.mid();

    if (c_f == 0.0) {
        return TaylorModel<T>(bound().inv());
    }

    // calculation using horner-scheme
    Vec<Iv> coeffs;
    coeffs.push_back(1.0/c_f);
    for (int i=1; i<=max_total_degree; i++) {
        coeffs.push_back(- coeffs.back() / c_f);
    }

    TaylorModel<T> ret = horner(f_bar,coeffs);

    // compute f_bar^(k+1)/c_f^(k+2) using interval arithmetics
    Iv tmp1 = pow(-1.0,max_total_degree+1)*f_bar.bound().power(max_total_degree+1)/c_f.power(max_total_degree+2);

    // compute 1.0 + \theta * f_bar / c_f using interval arithmetics
    Iv tmp2 = 1.0+Iv(0,1)*f_bar.bound()/c_f;

    ret.rest +=  tmp1/tmp2.power(max_total_degree+2);

    return ret;
}

/*!
 * The sine can be determined using the horner-scheme implemented in horner().
 * For more details, see \cite Makino2003.
 *
 * Note that in some cases the results differ depending on the value of
 * correct_trig. See correct_trig for more details.
 */
template <class T>
TaylorModel<T> TaylorModel<T>::sin() const
{
    if (max_total_degree<0)
        throw("sin(TaylorModel) - cannot be used with taylor-model of unlimited total degree!");

    // evaluate using interval arithmetics and intersect with [-1,1] (range of sine)
    Iv iv_ret = bound().sin().intersect(Iv(-1.0,1.0));
    if (iv_ret.isNaN())
        // if interval is NaN we still assume that the result is constrained to [-1,1]
        iv_ret = Iv(-1.0,1.0);

    // if pure interval-TaylorModel the result is again an interval
    if ( (poly.bound(domain) == 0.0) && !(rest==0.0) )
        return TaylorModel<T>(iv_ret);

    // if the interval remainder is unbounder of NaN the result is set to the interval [-1,1]
    if ( rest.isNaN() || (rest.lower() == -INFINITY) || (rest.upper() == INFINITY) )
        return TaylorModel<T>(iv_ret);

    // otherwise we proceede to calculate the result
    Iv c_f = cF().mid();
    TaylorModel<T> f_bar = *this - c_f.mid();

    // calculation using horner-scheme
    Vec<Iv> coeffs;
    double factorial=1;
    for (int i=0; i<=max_total_degree; i++) {
        if (i>1)
            factorial *= i;

        if (i%4==0)
            coeffs.push_back(c_f.sin()/factorial);
        else if (i%4==1)
            coeffs.push_back(c_f.cos()/factorial);
        else if (i%4==2)
            coeffs.push_back(-c_f.sin()/factorial);
        else
            coeffs.push_back(-c_f.cos()/factorial);
    }

    TaylorModel<T> ret = horner(f_bar,coeffs);

    // compute f_bar^(k+1)/(k+1)! using interval arithmetics
    factorial *= max_total_degree+1;
    Iv tmp1 = f_bar.bound().power(max_total_degree+1)/factorial;

    // compute J (from Makino, 2003)
    Iv tmp2 = c_f+Iv(0,1)*f_bar.bound();
    Iv tmp3;
    if ((max_total_degree+1)%4==0)
        tmp3 = tmp2.sin();
    else if ((max_total_degree+1)%4==1)
        tmp3 = tmp2.cos();
    else if ((max_total_degree+1)%4==2)
        tmp3 = -tmp2.sin();
    else
        tmp3 = -tmp2.cos();

    ret.rest += tmp1*tmp3;

    /* if correct_trig, then we replace the TaylorModel by an interval if the
     * remainder interval is NaN or larger than the complete interval dermined
     * using interval arithmetics.
     */
    if ( correct_trig && ( ret.rest.isNaN() || ( (iv_ret.lower() > ret.rest.lower()) && (iv_ret.upper() < ret.rest.upper() ) ) ) ) {
        return TaylorModel<T>(iv_ret);
    }

    return ret;
}

/*!
 * The cosine can be determined using the horner-scheme implemented in horner().
 * For more details, see \cite Makino2003.
 *
 * Note that in some cases the results differ depending on the value of
 * correct_trig. See correct_trig for more details.
 */
template <class T>
TaylorModel<T> TaylorModel<T>::cos() const
{
    if (max_total_degree<0)
        throw("cos(TaylorModel) - cannot be used with taylor-model of unlimited total degree!");

    // evaluate using interval arithmetics and intersect with [-1,1] (range of sine)
    Iv iv_ret = bound().cos().intersect(Iv(-1.0,1.0));
    if (iv_ret.isNaN())
        iv_ret = Iv(-1.0,1.0);

    // if pure interval-TaylorModel the result is again an interval
    if ( (poly.bound(domain) == 0.0) && !(rest==0.0) )
        return TaylorModel<T>(iv_ret);

    // if the interval remainder is unbounder of NaN the result is set to the interval [-1,1]
    if ( rest.isNaN() || (rest.lower() == -INFINITY) || (rest.upper() == INFINITY) ) {
        return TaylorModel<T>(iv_ret);
    }

    // otherwise we proceede to calculate the result
    Iv c_f = cF().mid();
    TaylorModel<T> f_bar = *this - c_f.mid();

    // calculation using horner-scheme
    Vec<Iv> coeffs;
    double factorial=1;
    for (int i=0; i<=max_total_degree; i++) {
        if (i>1)
            factorial *= i;

        if (i%4==3)
            coeffs.push_back(c_f.sin()/factorial);
        else if (i%4==0)
            coeffs.push_back(c_f.cos()/factorial);
        else if (i%4==1)
            coeffs.push_back(-c_f.sin()/factorial);
        else
            coeffs.push_back(-c_f.cos()/factorial);
    }

    TaylorModel<T> ret = horner(f_bar,coeffs);

    // compute f_bar^(k+1)/(k+1)! using interval arithmetics
    factorial *= max_total_degree+1;
    Iv tmp1 = f_bar.bound().power(max_total_degree+1)/factorial;

    // compute J (from Makino, 2003)
    Iv tmp2 = c_f+Iv(0,1)*f_bar.bound();
    Iv tmp3;
    if ((max_total_degree+1)%4==3)
        tmp3 = tmp2.sin();
    else if ((max_total_degree+1)%4==0)
        tmp3 = tmp2.cos();
    else if ((max_total_degree+1)%4==1)
        tmp3 = -tmp2.sin();
    else
        tmp3 = -tmp2.cos();

    ret.rest += tmp1*tmp3;

    /* if correct_trig, then we replace the TaylorModel by an interval if the
     * remainder interval is NaN or larger than the complete interval dermined
     * using interval arithmetics.
     */
    if ( correct_trig && ( ret.rest.isNaN() || ( (iv_ret.lower() > ret.rest.lower()) && (iv_ret.upper() < ret.rest.upper() ) ) ) ) {
        return TaylorModel<T>(iv_ret);
    }

    return ret;
}

/*!
 * The exponential function can be determined using the horner-scheme
 * implemented in horner(). For more details, see \cite Makino2003.
 */
template <class T>
TaylorModel<T> TaylorModel<T>::exp() const
{
    if (max_total_degree<0)
        throw("exp(TaylorModel) - cannot be used with taylor-model of unlimited total degree!");

    // if pure interval-TaylorModel the result is again an interval
    if ( (poly.bound(domain) == 0.0) && !(rest==0.0) )
        return TaylorModel<T>(rest.exp());

    Iv c_f = cF().mid();
    TaylorModel<T> f_bar = *this - c_f.mid();
    Iv c_f_exp(c_f.exp());

    // calculation using horner-scheme
    Vec<Iv> coeffs;
    double factorial=1;
    for (int i=0; i<=max_total_degree; i++) {
        if (i>1)
            factorial *= i;

        coeffs.push_back( (1.0/factorial) * c_f_exp);
    }

    TaylorModel<T> ret = horner(f_bar,coeffs);

    // compute f_bar^(k+1)/(k+1)! using interval arithmetics
    factorial *= max_total_degree+1;
    Iv tmp1 = f_bar.bound().power(max_total_degree+1)/factorial;

    // compute J (from Makino, 2003)
    Iv tmp2 = (Iv(0,1)*f_bar.bound()).exp();
    ret.rest += tmp1*tmp2;

    return ret;
}

/*!
 * Only implemented for non-negative integers. Naiively evaluated by multiplying
 * TaylorModel with itself.
 */
template <class T>
TaylorModel<T> TaylorModel<T>::power(size_t i) const
{
    if (i==0)
        return TaylorModel<T>(1.0);
    else if (i==1)
        return *this;
    else {
        TaylorModel ret = *this;
        for (size_t k=2; k<=i; k++)
            ret *= *this;
        return ret;
    }
}

template <class T>
bool TaylorModel<T>::operator==(double rhs) const
{
    return (bound() == rhs);
}

template <class T>
Iv TaylorModel<T>::cF() const
{
    return poly.constant();
}

/*!
 * See \cite Makino2003 for more details.
 */
template <class T>
TaylorModel<T> TaylorModel<T>::horner(const TaylorModel<T> &x, const Vec<Iv> &a) const
{
    TaylorModel<T> ret(0.0);

    for (int i=a.size()-1; i>=0; i--) {
        ret *= x;
        ret += a[i].mid();
        ret += a[i]-a[i].mid();
    }

    return ret;
}

/*!
 * Makes use of the Polynomial::cutoff().
 *
 * \todo Improve Polynomial::cutoff() to determine sharper bounds here.
 */
template <class T>
void TaylorModel<T>::cutoff()
{
    if (max_total_degree>=0)
        rest += poly.cutoff(max_total_degree,domain);
}

template class TaylorModel<double>;
template class TaylorModel<Iv>;

template <class T> std::complex<TaylorModel<T> > exp(const std::complex<TaylorModel<T> > &t) {

    return (t.real().exp() * std::complex<TaylorModel<T> >(t.imag().cos(),t.imag().sin()));
}

template std::complex<TaylorModel<double> > exp<double>(const std::complex<TaylorModel<double> > &t);
template std::complex<TaylorModel<Iv> > exp<Iv>(const std::complex<TaylorModel<Iv> > &t);

/*!
 * Explicitly considers the special cases where the imaginary part of the lhs or
 * rhs is zero for TaylorModels of type double. In our case this represents a
 * real TaylorModel, i.e., we assume that the result is zero, even when
 * multiplying this with NaN or infinity!
 */
template <> template <> std::complex<TaylorModel<double> > &std::complex<TaylorModel<double> >::operator*=(const std::complex<TaylorModel<double> > &rhs)
{
    if (imag()==0.0) {
        imag(real()*rhs.imag());
        real(real()*rhs.real());
    } else if (rhs.imag()==0.0) {
        imag(rhs.real()*imag());
        real(real()*rhs.real());
    } else {
        TaylorModel<double> r = real()*rhs.real() - imag()*rhs.imag();
        imag(real()*rhs.imag()+imag()*rhs.real());
        real(r);
    }
    return *this;
}

/*!
 * Explicitly considers the special cases where the imaginary part of the lhs or
 * rhs is zero for TaylorModels of type Iv. In our case this represents a real
 * TaylorModel, i.e., we assume that the result is zero, even when multiplying
 * this with NaN or infinity!
 *
 * \todo This special case should not be necessary anymore due to the fact that
 * the Iv-class can represent values that are "exactly_zero" in the sense of
 * this definition.
 */
template <> template <> std::complex<TaylorModel<Iv> > &std::complex<TaylorModel<Iv> >::operator*=(const std::complex<TaylorModel<Iv> > &rhs)
{
    if (imag()==0.0) {
        imag(real()*rhs.imag());
        real(real()*rhs.real());
    } else if (rhs.imag()==0.0) {
        imag(rhs.real()*imag());
        real(real()*rhs.real());
    } else {
        TaylorModel<Iv> r = real()*rhs.real() - imag()*rhs.imag();
        imag(real()*rhs.imag()+imag()*rhs.real());
        real(r);
    }
    return *this;
}

template std::complex<TaylorModel<double> > &std::complex<TaylorModel<double> >::operator*=(const std::complex<TaylorModel<double> > &rhs);
template std::complex<TaylorModel<Iv> > &std::complex<TaylorModel<Iv> >::operator*=(const std::complex<TaylorModel<Iv> > &rhs);

template <class T> std::complex<TaylorModel<T> > operator*(const std::complex<TaylorModel<T> > &t1,const std::complex<TaylorModel<T> > &t2)
{
    return (std::complex<TaylorModel<T> >(t1)*=t2);
}

template std::complex<TaylorModel<double> > operator*(const std::complex<TaylorModel<double> > &t1,const std::complex<TaylorModel<double> > &t2);
template std::complex<TaylorModel<Iv> > operator*(const std::complex<TaylorModel<Iv> > &t1,const std::complex<TaylorModel<Iv> > &t2);

template <class T>
std::ostream & operator<<(std::ostream &os, const TaylorModel<T>& t)
{
    std::ostringstream buf;

    buf << "Taylor-Model on domain " << t.domain << " with " << t.poly.val().size() << " coefficients:" << std::endl;
    buf << "poly        = " << t.poly << std::endl;
    buf << "poly_bound  = " << t.poly.bound(t.domain) << std::endl;
    buf << "rest        = " << t.rest << std::endl;
    buf << "bound       = " << t.bound() << std::endl;

    return (os << buf.str());
}

template std::ostream & operator<<(std::ostream &os, const TaylorModel<double>& t);
template std::ostream & operator<<(std::ostream &os, const TaylorModel<Iv>& t);
