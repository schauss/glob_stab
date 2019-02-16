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

#include "optimization.h"
#include "branch_and_bound.hpp"

#include "config.h"

/*!
 * Also initializes different statistics-counters and creates an
 * OptimizationResult Optimization::res. All of these are created as shared
 * pointers so they can be reused throughout the optimization/constraint
 * solving.
 */
Optimization::Optimization(size_t num_threads)
    : num_threads(num_threads),
      depth(0)
{
    res = boost::shared_ptr<OptimizationResult>(new OptimizationResult);

    id_count = boost::shared_ptr<size_t>(new size_t(0));
    id = (*id_count)++;

    split_B_count = boost::shared_ptr<size_t>(new size_t(0));
    split_Tm_count = boost::shared_ptr<size_t>(new size_t(0));
    split_par_count = boost::shared_ptr<size_t>(new size_t(0));

    fraction_done = boost::shared_ptr<double>(new double(0.0));
}

void Optimization::optimize(const Function &cost)
{
    fcn = boost::shared_ptr<Function>(new Function(cost));
    (*fraction_done) = 1.0 - fcn->parameters().relativeVolume();

    BranchAndBound bab(num_threads);
    bab.run(*this,&Optimization::doOptimize);
}

void Optimization::constrain(const Function &constraint)
{
    fcn = boost::shared_ptr<Function>(new Function(constraint));
    (*fraction_done) = 1.0 - fcn->parameters().relativeVolume();

    BranchAndBound bab(num_threads);
    bab.run(*this,&Optimization::doConstrain);
}

std::vector<Optimization> Optimization::doOptimize()
{
    return std::vector<Optimization>();
}

std::vector<Optimization> Optimization::doConstrain()
{
    if (res->isDone(fcn->parameters())) {
        (*fraction_done) += fcn->parameters().relativeVolume();
        printProgress("already done ");
        return std::vector<Optimization>();
    }

    if (!fcn->isEvaluated()) {
        fcn->evaluate();
    }

    if (fcn->excludesZero()) {

        (*fraction_done) += fcn->parameters().relativeVolume();

        if ( !(fcn->doValueSet()) && (fcn->value()<=0.0).any() ) {
            /* TODO
             * if value <= 0 and !(value <0) this means that upper bound of value == 0
             * => actually we would have to check if upper bound of value is sharp, if it is
             *    the system is marginally stable, if it isn't we cannot be sure
             * => in practice an upper bound of exactly zero for a_0/a_1 with s=0
             *    is unlikely if it is not sharp!
             */
            if ( !(fcn->hasMarginal()) && !((fcn->value()<0.0).any()) ) {
                res->addMarginal(fcn->parameters());
                printProgress("MARGINALLY STABLE!");
            } else {
                res->setDone(fcn->parameters());
                res->addViolated(fcn->parameters());
                printProgress("excludes+done");
            }
        } else {
            printProgress("excludes zero");
        }

        return std::vector<Optimization>();
    } else if (fcn->includesZero()) {

        if (fcn->parameters().minFractionReached()) {
            res->setDone(fcn->parameters());
            if (fcn->doValueSet() && !(fcn->doValueSetBoundary()) )
                res->addViolated(fcn->parameters());
            else
                res->addBoundary(fcn->parameters());

            (*fraction_done) += fcn->parameters().relativeVolume();
            printProgress("includes zero");
            return std::vector<Optimization>();
        } else {
            printProgress("partial      ");
            return splitPar();
        }
    } else if (fcn->partiallyIncludesZero1()) {
        printProgress("part inc one ");
        return splitGradB();
    } else if (fcn->partiallyIncludesZero2()) {

        if (fcn->isEvaluatedFromEquation()) {
            printProgress("part inc-spl ");
            return splitGradTm();
        } else {
            printProgress("part inc-reev");
            fcn->clear();
            return doConstrain();
        }
    } else {
        if (!fcn->restIsFinite()) {
            printProgress("unknown NAN  ");
            return splitGradTm();
        } else {
            printProgress("unknown      ");
            return splitGradB();
        }
    }
}

std::vector<Optimization> Optimization::splitGradB() const
{
    (*split_B_count)++;

    int split_dir = 0;

    Vec<Vec<double> > dx = fcn->dpoly_dparam();
    Vec<double> dx_max;
    Vec<size_t> dir_max;
    for(size_t i=0; i< dx.size(); i++) {
        double tmp = dx[i].maxVal(&split_dir);
        if ( (tmp <= 0.0) || isnan(tmp))
            continue;

        dx_max.push_back(tmp);
        dir_max.push_back(split_dir);
    }
    if (dx_max.empty()) {
        // if bernstein derivative is zero, evaluate tm-derivative
        return splitGradTm();
    }

    dx_max.maxVal(&split_dir);
    if (CONFIG_VAR(verbosity,int) > 0)
        std::cout << "gradient: " << dx << " -> split dir: " << dir_max[split_dir] << std::endl;
    return splitAt(dir_max[split_dir]);
}

std::vector<Optimization> Optimization::splitGradTm() const
{
    (*split_Tm_count)++;
    int split_dir = 0;

    Vec<Vec<double> > dx = fcn->drest_dparam();
    Vec<double> dx_min;
    Vec<size_t> dir_min;
    for(size_t i=0; i< dx.size(); i++) {

        double tmp = INFINITY;
        for (size_t j=0; j<dx[i].size(); j++) {
            if ( (fcn->parameters().value()[j].width() == 0.0) || ( (fcn->parameters().name(j) == "sigma") && fcn->doValueSetBoundary() ) )
                continue;
            else if ( (dx[i][j] < tmp) || ( (dx[i][j] == tmp) && (fcn->parameters().fraction(j) > fcn->parameters().fraction(split_dir)) ) ) {
                tmp = dx[i][j];
                split_dir = j;
            }
        }
        //double tmp = dx[i].minVal(&split_dir);

        //if ( (tmp <= 0.0) || isnan(tmp))
        //    continue;

        dx_min.push_back(tmp);
        dir_min.push_back(split_dir);
    }
    if (dx_min.empty()) {
        std::cout << "function:" << std::endl << *fcn << std::endl;
        std::cout << "gradient: " << dx << std::endl;
        throw("Optimization::splitGradTm - empty derivative!");
    }

    dx_min.maxVal(&split_dir);
    if (CONFIG_VAR(verbosity,int) > 0)
        std::cout << "gradient: " << dx << " -> split dir: " << dir_min[split_dir] << std::endl;

    fcn->clear();
    return splitAt(dir_min[split_dir]);
}

std::vector<Optimization> Optimization::splitPar() const
{
    (*split_par_count)++;

    return splitAt(fcn->parameters().maxRelFractionDir());
}

/*!
 * \todo Currently this function is always called without specifying mid_point,
 * i.e., mid_point=0.5 is used! This could lead to non-convergence, e.g., if a
 * root touches the boundary at 1/2, 1/4, 1/8, ... of the parameter-interval!
 * * Could be solved by randomly varying the split-point, say between 0.4 - 0.6
 * * Does this have any disadvantage? We could sometimes reach q_min one depth
 *   later, but that shouldn't cause a large overhead!
 * * Does not solve the problem when the root touches the boundary on the edge
 *   of the original interval, i.e., at 0 or 1.
 */
std::vector<Optimization> Optimization::splitAt(size_t direction, double mid_point) const
{
    std::vector<boost::shared_ptr<Function> > f = fcn->splitAt(direction,mid_point);

    std::vector<Optimization> ret(f.size(),*this);
    for (size_t i=0; i<f.size(); i++) {
        ret[i].fcn = f[i];
        ret[i].depth = depth+1;
        ret[i].id = (*id_count)++;
    }

    return ret;
}

void Optimization::printProgress(const std::string &state) const
{
    std::ostringstream buf;
    buf << "ID: " << std::setw(5) << id << ", depth " << std::setw(4) << depth << " - " << state  << "  " << std::fixed << std::setprecision(3) << (*fraction_done)*100.0 << " \% done" << std::endl << std::flush;
    if (CONFIG_VAR(verbosity,int) > 0)
        buf << "function:" << std::endl << *fcn << std::endl;
    std::cout << buf.str();
}


std::ostream & operator<<(std::ostream &os, const Optimization& opt)
{
    os << "Count: " << *(opt.id_count) << std::endl;
    os << "Splits: " << std::endl;
    os << " - Bernstein: " << *(opt.split_B_count) << std::endl;
    os << " - Taylor:    " << *(opt.split_Tm_count) << std::endl;
    os << " - Parameter: " << *(opt.split_par_count) << std::endl;
    return os;
}
