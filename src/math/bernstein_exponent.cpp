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

#include "bernstein_exponent.h"

#include <boost/math/special_functions/binomial.hpp>
#include <cstdlib>

std::map<std::vector<uint8_t>,BernsteinExponent*> BernsteinExponent::instance;
boost::mutex  BernsteinExponent::instance_mutex;

BernsteinExponent *BernsteinExponent::getInstance(const std::vector<uint8_t> &exp)
{
    BernsteinExponent *ex;

    instance_mutex.lock();

    std::map<std::vector<uint8_t>,BernsteinExponent*>::const_iterator it = instance.find(exp);
    if (it == instance.end()) {
        ex = new BernsteinExponent(exp);
        instance[exp] = ex;
    } else {
        ex = it->second;
    }

    instance_mutex.unlock();
    return ex;
}

void BernsteinExponent::clear()
{
    instance_mutex.lock();
    std::map<std::vector<uint8_t>,BernsteinExponent*>::const_iterator it;
    for (it = instance.begin(); it!= instance.end(); it++)
        delete it->second;
    instance.clear();
    instance_mutex.unlock();
}

BernsteinExponent::BernsteinExponent(const std::vector<uint8_t> &exp)
    : max_degree(exp)
{
    std::cout << "Calculating Exponent for max_degree " << Vec<int>(max_degree) << std::endl;

    /*
       calculate steps from one element of this dimension to next
       --> cumulative product of [1 max_degree+1]
    */
    step.push_back(1);
    Vec<size_t> tmp = max_degree;
    tmp += 1;
    step.push_back(tmp.cumProd());

    // step is one element to large... last element of step is number of ids
    max_id = step.back();
    step.pop_back();

    //std::cout << "steps = " << step << std::endl;
    //std::cout << "max_id = " << max_id << std::endl;

    /* temporary variables for following computations... */
    uint8_t num_var = exp.size();
    uint8_t order = max_degree.maxVal();

    vec_id = Vec<uint8_t>(max_id*num_var,0);
    vertex_id = Vec<uint8_t>(pow(2.0,num_var)*boost::math::binomial_coefficient<double>(num_var,0),0);
    lin_id = Vec<uint8_t>(num_var+1,0);
    uint8_t vertex_id_count = 0; // extract the id's for vertexes as I currently can't come up with the (simple) formula to compute them!
    uint8_t lin_id_count = 0;
    bool is_vertex;
    uint16_t dim_sum_count;

    ldiv_t divresult;
    for (unsigned int id = 0; id<max_id; id++) {
        divresult.rem = (long)id;
        dim_sum_count = 0;
        is_vertex = true;
        for (int i = (int)num_var-1; i>=0; i--) {
            divresult = ldiv(divresult.rem, (long)step[i]);
            vec_id[id*num_var+i] = divresult.quot;
            dim_sum_count += divresult.quot;
            if ( (divresult.quot != max_degree[i]) && (divresult.quot != 0) )
                is_vertex = false;
        }
        // extract id's of vertexes
        if ( is_vertex ) {
            //std::cout << "Vertex with id " << id << "-> " << Vec<int>(toVec(id)) << std::endl;
            vertex_id[vertex_id_count++] = id;
        }

        // extract id's of linear part
        if (dim_sum_count <= 1)
            lin_id[lin_id_count++] = id     ;
    }

    //printf("ID-vectors set up!\n");

    Vec<double> binom((int)(order+1)*(int)(order+1),0.0);
    // calculate binomial coefficients
    for (int i=0; i<=order; i++) {
        for (int j=0; j<=order; j++) {
            if (j <= i)
                binom[i*(order+1)+j] = boost::math::binomial_coefficient<double>(i,j);
        }
    }

    //printf("Binomial coefficients calculated!\n");

    uint8_t N[num_var];
    for(size_t i=0; i<num_var; i++)
        N[i] = max_degree[i];

    //printf("Coefficient N initialized!\n");

    binom_vec = Vec<double>(max_id,0.0);
    // calculate vector if 1/binomialCoefficients
    for (unsigned int i = 0; i<max_id; i++) {
        const uint8_t *J = toArray(i);
        binom_vec[i] = 1.0 / binomCoeffVec(N,J,binom,order);
    }
    //printf("Binomial coefficient vector calculated!\n\n");
}

BernsteinExponent::~BernsteinExponent()
{
}

const uint8_t *BernsteinExponent::toArray(unsigned int id) const
{
    return (vec_id.data()+id*numVar());
}

Vec<uint8_t> BernsteinExponent::toVec(unsigned int id) const
{
    Vec<uint8_t> v;
    for (size_t i=0; i<numVar(); i++)
        v.push_back(vec_id[id*numVar()+i]);

    return v;
}

unsigned int BernsteinExponent::toId(const uint8_t exponents[]) const
{
    unsigned int id = 0;
    for (size_t i=0; i<numVar(); i++) {
        id += exponents[i]*step[i];
    }

    return id;
}

unsigned int BernsteinExponent::toId(const Vec<uint8_t> &exponents) const
{
    if (exponents.size() != numVar())
        throw("Exponent::toId - size missmatch!");

    unsigned int id = 0;
    for (size_t i=0; i<numVar(); i++) {
        id += exponents[i]*step[i];
    }

    return id;
}

unsigned int BernsteinExponent::maxId() const
{
    return max_id;
}

bool BernsteinExponent::isVertex(unsigned int id) const
{
    for (unsigned int i = 0; i < vertex_id.size(); i++) {
        if (vertex_id[i] == id)
            return true;
    }
    return false;
}

// only still used in initialization!
double BernsteinExponent::binomCoeffVec(const uint8_t n[], const uint8_t k[], const Vec<double> &binom, uint8_t order)
{
    double res = 1.0;

    for (size_t i=0; i<numVar(); i++)
        res *= binom[n[i]*(order+1)+k[i]];

    return res;
}
