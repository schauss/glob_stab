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

#include "optimization_result.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>

OptimizationResult::OptimizationResult()
{
    Vec<size_t> end_idx = Parameters().endIdx();
    done = ndArray<bool>(end_idx);
}

OptimizationResult::~OptimizationResult()
{

}

void OptimizationResult::addBoundary(const Parameters &param)
{
    addResult(param,constraint_boundary);
}

void OptimizationResult::addSatisfied(const Parameters &param)
{
    addResult(param,constraint_satisfied);
}

void OptimizationResult::addViolated(const Parameters &param)
{
    addResult(param,constraint_violated);
}

void OptimizationResult::addMarginal(const Parameters &param)
{
    addResult(param,constraint_marginal);
}

void OptimizationResult::addResult(const Parameters &param, Vec<Parameters> &res)
{
    typename std::map<std::vector<int>,std::vector<std::vector<int> > >::const_iterator it = final_test.find(Vec<int>(param.startIdx()));
    boost::mutex::scoped_lock l(store_mutex);
    if (it == final_test.end())
        res.push_back(param.indexExtents());
    else {
        for (size_t i=0; i< it->second.size(); i++) {
            res.push_back(Parameters(Vec<size_t>(it->second[i]), Vec<size_t>(it->second[i])+1));
        }
    }
}

bool OptimizationResult::hasMarginal() const
{
    return (!constraint_marginal.empty());
}

bool OptimizationResult::storeAsy(std::string filename) const
{
    std::ofstream hnd(filename.c_str());
    if (!hnd.good()) {
        throw("OptimizationResult::storeAsy - Error opening output file!");

    }

    hnd << std::setprecision(15);

    Parameters par;

    // find parameters which are sampled
    Vec<size_t> idx;
    for (size_t i = 0; i<par.size(); i++) {
        if (par.minFraction(i) != 0)
            idx.push_back(i);
    }

    if (idx.size() != 2) {
        std::cout << "Cannot store asy for " << idx.size() << " sampled parameters (only 2)" << std::endl;
        hnd.close();
        return false;
    }

    // write header
    hnd << "import graph;" << std::endl;
    hnd << "include glob_stab_fig;" << std::endl << std::endl;

    Iv a,b;
    //store complete region as stable, this must then be overlayd with unstable region!
    a = par.value()[idx[0]];
    b = par.value()[idx[1]];
    hnd << "filldraw(box((" << a.lower() << "," << b.lower() << "),(" << a.upper() << "," << b.upper() << ")),white,white+0);" << std::endl;

    for (size_t i = 0; i<constraint_satisfied.size(); i++) {
        a = constraint_satisfied[i].value()[idx[0]];
        b = constraint_satisfied[i].value()[idx[1]];
        hnd << "filldraw(box((" << a.lower() << "," << b.lower() << "),(" << a.upper() << "," << b.upper() << ")),green,gray+0);" << std::endl;
    }

    for (size_t i = 0; i<constraint_marginal.size(); i++) {
        a = constraint_marginal[i].value()[idx[0]];
        b = constraint_marginal[i].value()[idx[1]];
        hnd << "filldraw(box((" << a.lower() << "," << b.lower() << "),(" << a.upper() << "," << b.upper() << ")),lightgray,lightgray+0);" << std::endl;
    }

    for (size_t i = 0; i<constraint_boundary.size(); i++) {
        a = constraint_boundary[i].value()[idx[0]];
        b = constraint_boundary[i].value()[idx[1]];
        hnd << "filldraw(box((" << a.lower() << "," << b.lower() << "),(" << a.upper() << "," << b.upper() << ")),black,black+0);" << std::endl;
    }

    for (size_t i = 0; i<constraint_violated.size(); i++) {
        a = constraint_violated[i].value()[idx[0]];
        b = constraint_violated[i].value()[idx[1]];
        hnd << "filldraw(box((" << a.lower() << "," << b.lower() << "),(" << a.upper() << "," << b.upper() << ")),mediumgray,mediumgray+0);" << std::endl;
    }

    hnd << "xaxis(\"$" << par.name(idx[0]) << "$\",Bottom,LeftTicks,above=true);" << std::endl;
    hnd << "yaxis(\"$" << par.name(idx[1]) << "$\",Left,RightTicks,above=true);" << std::endl;

    hnd.close();

    return true;
}

void OptimizationResult::storeMatlab(std::string filename) const
{
    Parameters par;

    // find parameters which are sampled
    Vec<size_t> idx;
    for (size_t i = 0; i<par.size(); i++) {
        if (par.minFraction(i) != 0)
            idx.push_back(i);
    }

    if (idx.empty()) {
        std::cout << "No sampled parameters! -> not saving file!" << std::endl;
        return;
    }

    std::ofstream hnd(filename.c_str());
    if (!hnd.good()) {
        throw("Transferfunction::store - Error opening output file!");

    }

    hnd << std::setprecision(15);

    std::cout << "Saving results to file " << filename << std::endl;

    hnd << "parameters = {";
    for (size_t i = 0; i<idx.size(); i++) {
        hnd << "'" << par.name(idx[i]) << "' ";
    }
    hnd << "};";

    hnd << std::endl << std::endl;

    hnd << "ranges = zeros(" << idx.size() << ", 2);" << std::endl;
    for (size_t j = 0; j<idx.size(); j++) {
        Iv val = par.value()[idx[j]];
        hnd << "ranges(" << j+1 << ", 1) = " << val.lower() << ";" << std::endl;
        hnd << "ranges(" << j+1 << ", 2) = " << val.upper() << ";" << std::endl;
    }
    hnd << std::endl;

    if (!constraint_satisfied.empty()) {
        hnd << "constraint_satisfied = zeros(" << constraint_satisfied.size() << ", " << idx.size() << ", 2);" << std::endl;
        for (size_t i = 0; i<constraint_satisfied.size(); i++) {
            for (size_t j = 0; j<idx.size(); j++) {
                Iv val = constraint_satisfied[i].value()[idx[j]];
                hnd << "constraint_satisfied(" << i+1 << ", " << j+1 << ", 1) = " << val.lower() << ";" << std::endl;
                hnd << "constraint_satisfied(" << i+1 << ", " << j+1 << ", 2) = " << val.upper() << ";" << std::endl;
            }
        }
        hnd << std::endl;
    }

    if (!constraint_marginal.empty()) {
        hnd << "constraint_marginal = zeros(" << constraint_marginal.size() << ", " << idx.size() << ", 2);" << std::endl;
        for (size_t i = 0; i<constraint_marginal.size(); i++) {
            for (size_t j = 0; j<idx.size(); j++) {
                Iv val = constraint_marginal[i].value()[idx[j]];
                hnd << "constraint_marginal(" << i+1 << ", " << j+1 << ", 1) = " << val.lower() << ";" << std::endl;
                hnd << "constraint_marginal(" << i+1 << ", " << j+1 << ", 2) = " << val.upper() << ";" << std::endl;
            }
        }
        hnd << std::endl;
    }

    if (!constraint_boundary.empty()) {
        hnd << "constraint_boundary = zeros(" << constraint_boundary.size() << ", " << idx.size() << ", 2);" << std::endl;
        for (size_t i = 0; i<constraint_boundary.size(); i++) {
            for (size_t j = 0; j<idx.size(); j++) {
                Iv val = constraint_boundary[i].value()[idx[j]];
                hnd << "constraint_boundary(" << i+1 << ", " << j+1 << ", 1) = " << val.lower() << ";" << std::endl;
                hnd << "constraint_boundary(" << i+1 << ", " << j+1 << ", 2) = " << val.upper() << ";" << std::endl;
            }
        }
        hnd << std::endl;
    }

    if (!constraint_violated.empty()) {
        hnd << "constraint_violated = zeros(" << constraint_violated.size() << ", " << idx.size() << ", 2);" << std::endl;
        for (size_t i = 0; i<constraint_violated.size(); i++) {
            for (size_t j = 0; j<idx.size(); j++) {
                Iv val = constraint_violated[i].value()[idx[j]];
                hnd << "constraint_violated(" << i+1 << ", " << j+1 << ", 1) = " << val.lower() << ";" << std::endl;
                hnd << "constraint_violated(" << i+1 << ", " << j+1 << ", 2) = " << val.upper() << ";" << std::endl;
            }
        }
        hnd << std::endl;
    }

    hnd.close();
}

bool OptimizationResult::isDone(const Parameters &param) const
{
    return (done(param.startIdx(),param.endIdx()) == true);
}

bool OptimizationResult::anyDone(const Parameters &param) const
{
    bool ret = done(param.startIdx(),param.endIdx()).any(true);
/*
    if (ret) {
        std::cout << "\nIDX-Range: ";
        for(size_t i=0; i<start_idx.size();i++)
            std::cout << start_idx[i] << "-" << end_idx[i]-1 << " ";
        std::cout << std::endl;

        std::cout << "Done:\n";
        std::vector<std::vector<size_t> > d_vec = done(start_idx,end_idx).find(true);
        for(size_t i=0; i<d_vec.size();i++) {
            for(size_t j=0; j<d_vec[i].size();j++)
                std::cout << d_vec[i][j] << " ";
            std::cout << std::endl;
        }
    }
*/
    return ret;
}

void OptimizationResult::setDone(const Parameters &param)
{
    done(param.startIdx(),param.endIdx()) = true;
}

Parameters OptimizationResult::shrinkNotDone(const Parameters &param) const
{
    if (isDone(param))
        throw("OptimizationResult::shrinkNotDone() - isDone!\n");

    std::vector<size_t> start,end;
    bool found_not_done = done(param.startIdx(),param.endIdx()).boundingBox(false,start,end);

    if (!found_not_done)
        throw("OptimizationResult::shrinkNotDone() - all done!\n");

    Parameters bounding_box(start,end);
    return param.intersect(bounding_box);
}

Vec<Parameters> OptimizationResult::finalTestParameters()
{
    Parameters original;
    std::vector<std::vector<size_t> > not_done_vec = done(original.startIdx(),original.endIdx()).find(false);
    std::vector<std::vector<size_t> > done_vec = done(original.startIdx(),original.endIdx()).find(true);

    std::set<std::vector<int> > missing;

    for (size_t i=0; i< not_done_vec.size(); i++) {
        missing.insert(Vec<int>(not_done_vec[i]));
    }

    while (!missing.empty()) {
        std::vector<std::vector<int> > region(1,*(missing.begin()));
        missing.erase(missing.begin());
        addMissing(region,missing);
        final_test[findFurthest(region,done_vec)] = region;
    }

    std::cout << std::endl << "Checking " << final_test.size() << " disjoint regions!" << std::endl;
    Vec<Parameters> ret;
    for (typename std::map<std::vector<int>,std::vector<std::vector<int> > >::const_iterator it = final_test.begin(); it!=final_test.end(); it++) {
        Parameters mid(Vec<size_t>(it->first), Vec<size_t>(it->first)+1);
        mid.scale(0.0);
        ret.push_back(mid);
    }
    return ret;
}

void OptimizationResult::addMissing(std::vector<std::vector<int> > &region, std::set<std::vector<int> > &missing) const
{
    typename std::set<std::vector<int> >::iterator it;

    for (size_t idx=0; idx<region.size(); idx++) {
        std::vector<int> mid(region[idx]);
        for (size_t i=0; i<mid.size(); i++) {
            std::vector<int> tmp(mid);

            tmp[i]++;
            it = missing.find(tmp);
            if ( it != missing.end()) {
                region.push_back(*it);
                missing.erase(it);
            }

            tmp[i] -= 2;
            it = missing.find(tmp);
            if ( it != missing.end()) {
                region.push_back(*it);
                missing.erase(it);
            }
        }
    }
}

std::vector<int> OptimizationResult::findFurthest(const std::vector<std::vector<int> > &region, const std::vector<std::vector<size_t> > &done_vec) const
{
    // find maximum and minimum index for each direction
    size_t min_dist = 0;
    size_t min_idx = 0;

    Vec<size_t> idx;

    Parameters p;
    for (size_t i=0; i<p.size(); i++) {
        if ( (p.minFraction(i) != 0.0) && (p.minFraction(i) != 1.0) )
            idx.push_back(i);
    }

    for (size_t i=0; i<region.size(); i++) {
        size_t min_dist_i = 0;
        for (size_t j=0; j<done_vec.size(); j++) {
            size_t dist_j = 0;
            for (size_t k=0; k<idx.size(); k++) {
                int delta = region[i][idx[k]]-done_vec[j][idx[k]];
                dist_j += delta*delta;
            }
            if ( (j==0) || (dist_j < min_dist_i) ) {
                min_dist_i = dist_j;
            }
        }
        if ( (i==0) || (min_dist_i > min_dist) ) {
            min_dist = min_dist_i;
            min_idx = i;
        }
    }

    std::cout << "Checking region " << Vec<int>(region[min_idx]) << " with minimum distance: " << min_dist << std::endl;
    return region[min_idx];
}

Vec<Parameters> OptimizationResult::getAllNotDone() const
{
    return getNotDone(Parameters());
}

Vec<Parameters> OptimizationResult::getNotDone(const Parameters &param) const
{
    std::vector<std::vector<size_t> > not_done_vec = done(param.startIdx(),param.endIdx()).find(false);

    Vec <Parameters> ret;
    for (size_t i=0; i<not_done_vec.size(); i++) {
        ret.push_back(Parameters(not_done_vec[i],Vec<size_t>(not_done_vec[i])+1));
    }
    return ret;
}

std::ostream & operator<<(std::ostream &os, const OptimizationResult& res)
{
    if (res.constraint_boundary.empty() && res.constraint_violated.empty()) {
        os << "Transfer function is stable for complete parameter range!\n";
        return os;
    }

    os << "Transfer function is NOT stable for complete parameter range...\n";

    os << "\n---------------------------------------------------------------------------\n";
    os << " stable for\n";
    os << "---------------------------------------------------------------------------\n\n";
    for (size_t i = 0; i<res.constraint_satisfied.size(); i++) {
        os << res.constraint_satisfied[i].print();
    }

    os << "\n---------------------------------------------------------------------------\n";
    os << " marginally stable for\n";
    os << "---------------------------------------------------------------------------\n\n";
    for (size_t i = 0; i<res.constraint_marginal.size(); i++) {
        os << res.constraint_marginal[i].print();
    }
    os << "\n---------------------------------------------------------------------------\n";


    os << "\n---------------------------------------------------------------------------\n";
    os << " on boundary for\n";
    os << "---------------------------------------------------------------------------\n\n";
    for (size_t i = 0; i<res.constraint_boundary.size(); i++) {
        os << res.constraint_boundary[i].print();
    }
    os << "\n---------------------------------------------------------------------------\n";

    os << "\n---------------------------------------------------------------------------\n";
    os << " unstable for\n";
    os << "---------------------------------------------------------------------------\n\n";
    for (size_t i = 0; i<res.constraint_violated.size(); i++) {
        os << res.constraint_violated[i].print();
    }
    os << "\n---------------------------------------------------------------------------\n";

    return os;
}
