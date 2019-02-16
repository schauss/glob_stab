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

#ifndef OPTIMIZATION_RESULT_H
#define OPTIMIZATION_RESULT_H

#include <boost/thread/mutex.hpp>
#include <map>
#include <string>

#include "parameters.h"
#include "vec.hpp"

/*!
 * \brief The class OptimizationResult stores optimization results obtained in
 * the Optimization class.
 *
 * In its current form it actually stores results of a constraint solving
 * problem in the parameter space.
 *
 * The functions addBoundatay(), addSatisfied(), and addViolated() are used to
 * set whether, for the given parameter set, we are on a boundary, the
 * constraints are satisfied, or violated. In addition, a special case only
 * needed for stability analysis is considered via the addMarginal() function
 * which allows to specify that the parameter set is marginally stable.
 *
 * This class is also used to store regions which do not have to be considered
 * anymore and check this information using the functions setDone(), isDone(),
 * and anyDone().
 *
 * The resulting stable regions can either stored using storeMatlab() (MATLAB
 * m-file) or using storeAsy() (Asymptote plot). In addition, a human-readable
 * stream can be output using the stream operator
 * std::ostream & operator<<(std::ostream &os, const OptimizationResult& res).
 */
class OptimizationResult
{
public:
    /*!
     * \brief Construct a new instance and initialize OptimizationResult::done.
     * Call after initial Parameters have been set!
     *
     * The number of dimensions and length of each dimension is set according to
     * the value in Paramters().end_idx(), i.e., this depends on the number of
     * interval paramters in Parameters and on the set resolution for each
     * paramter.
     */
    OptimizationResult();
    //! Destructor.
    ~OptimizationResult();

    /*!
     * \brief Mark \a param as boundary, i.e., add \a param to
     * OptimizationResult::constraint_boundary.
     *
     * Calls addResult() to do this, see that function for more information.
     */
    void addBoundary (const Parameters &param);
    /*!
     * \brief Mark \a param as satisfying constraint, i.e., add \a param to
     * OptimizationResult::constraint_satisfied.
     *
     * Calls addResult() to do this, see that function for more information.
     */
    void addSatisfied(const Parameters &param);
    /*!
     * \brief Mark \a param as violating constraint, i.e., add \a param to
     * OptimizationResult::constraint_violated.
     *
     * Calls addResult() to do this, see that function for more information.
     */
    void addViolated (const Parameters &param);
    /*!
     * \brief Mark \a param as marginally stable, i.e., add \a param to
     * OptimizationResult::constraint_marginal.
     *
     * Calls addResult() to do this, see that function for more information.
     */
    void addMarginal (const Parameters &param);
    /*!
     * \brief Add result to one of the result vectors.
     * \param param Parameters to add to result vectors.
     * \param res One of
     * * OptimizationResult::constraint_boundary
     * * OptimizationResult::constraint_satisfied
     * * OptimizationResult::constraint_violated
     * * OptimizationResult::constraint_marginal
     *
     * When the set of disjoint regions OptimizationResult::final_test is empty
     * then \a param, enlarged to the desired resolution for boundary mapping
     * by calling Paramters::indexEntents(), is appeded to \a res.
     * If OptimizationResult::final_test is not empty then we will certainly
     * find param.startIndex() in OptimizationResult::final_test. In this case
     * we add all Parameter sets in this disjoint region to \a res.
     *
     * Uses OptimizationResult::store_mutex for thread-safety.
     *
     * \todo No need for this to be public!
     */
    void addResult(const Parameters &param, Vec<Parameters> &res);

    /*!
     * \brief Is the system marginally stable for some paramter?
     *
     * Checks whether OptimizationResult::constraint_marginal is not empty,
     * i.e., has addMarginal() been called?
     *
     * \todo Limited to the special case of stability analysis! Doesn't really
     * fit into this general purpose class for constraint-solving and
     * optimization. The same holds for addMarginal().
     */
    bool hasMarginal() const;

    /*!
     * \brief Store results as 2d [Asymptote](http://asymptote.sourceforge.net)
     * plot with given \a filename.
     * \return False, if the number of rastered parameters is not equal to two.
     * True otherwise, i.e., if the result was stored.
     *
     * Only implemented for 2d-plots, i.e., the desired resolution must be
     * specified for exactly two paramters!
     */
    bool storeAsy(std::string filename) const;
    /*!
     * \brief Store n-D results as MATLAB m-file with given \a filename.
     *
     * The different constraints are stored as 3d-array of dimension (m, n, 2)
     * where \em m is the number of regions for this constraint and \em n is the
     * number of rastered parameters and the lower and upper bound is stored for
     * each interval.
     */
    void storeMatlab(std::string filename) const;

    /*!
     * \brief Check whether \a param is completely done.
     *
     * Checks whether all values in the range between \a param.startIndex() and
     * param.endIndex() are set to true in OptimizationResult::done.
     */
    bool isDone(const Parameters &param) const;
    /*!
     * \brief Check whether any part of \a param is done.
     *
     * Checks whether any values in the range between \a param.startIndex() and
     * param.endIndex() are set to true in OptimizationResult::done.
     */
    bool anyDone(const Parameters &param) const;
    /*!
     * \brief Mark the region \a param as done.
     *
     * Set all values in the range between \a param.startIndex() and
     * param.endIndex() to true in OptimizationResult::done.
     */
    void setDone(const Parameters &param);

    //! Determine bounding box in parameter set \a param which is not done
    Parameters shrinkNotDone(const Parameters &param) const;
    /*!
     * \brief Get vector of all not-done regions in the initial parameter set
     * using getNotDone().
     *
     * \todo Could simply use Parameters() as default argument in getNotDone()
     * and remove this function.
     */
    Vec<Parameters> getAllNotDone() const;
    /*!
     * \brief Get vector of all not-done regions in the parameter set \a param.
     *
     * Note, that these regions have the same size as the desired resolution,
     * i.e., no aggregation of connected regions to bigger hyper-boxes is
     * implemented.
     */
    Vec<Parameters> getNotDone(const Parameters &param) const;

    /*!
     * \brief Determine disjoint regions separated by boundaries.
     * \return Vector of one parameterization per disjoint region.
     *
     * Finds disjoint regions in the parameter set using addMissing(). For each
     * of these disjoint regions the parameter set with the largest distance to
     * the boundaries is then determined using findFurthest() and this is used
     * as key in the map OptimizationResult::final_test where the disjoint
     * regions are stored (which in turn is used in addResult() to determine
     * which regions are affected by a result).
     */
    Vec<Parameters> finalTestParameters();
    /*!
     * \brief Finds connected regions, comparable to the *magic wand* in
     * graphics editors but in n dimensions.
     *
     * For all indices given in \a region find all indices in \a missing which
     * are directly connected to these start indices. Add them to \a region and
     * remove them from \a missing. Do this until we cannot find anymore points
     * in \a missing which are connected to points in \a region.
     *
     * \todo Rename to something understandable.
     */
    void addMissing(std::vector<std::vector<int> > &region, std::set<std::vector<int> > &missing) const;
    /*!
     * \brief Find index in \a region with largest distance from all entries in
     * \a done_vec.
     *
     * Currently this can return regions on the edges or faces of the original
     * parameter set as we only evaluate the distance from the done_vec and
     * not from the minimum and maximum index. This may not be ideal in some
     * cases as this could be close to a boundary we can simply not see because
     * it is outside our region of interest.
     *
     * \todo We could easily also take the minimum and maximum index into
     * account here. This could be selected using a parameter.
     */
    std::vector<int> findFurthest(const std::vector<std::vector<int> > &region, const std::vector<std::vector<size_t> > &done_vec) const;

private:
    //! Stores parameter sets which are on the boundary.
    Vec<Parameters> constraint_boundary;
    //! Stores parameter sets for which the constraints are satisfied.
    Vec<Parameters> constraint_satisfied;
    //! Stores parameter sets for which the constraints are violated.
    Vec<Parameters> constraint_violated;
    //! Stores parameter sets which are marginally stable.
    Vec<Parameters> constraint_marginal;

    //! Map of disjoint regions (from test-parameter to whole region)
    std::map<std::vector<int>,std::vector<std::vector<int> > > final_test;

    //! Regions which are done and don't have to be evaluated anymore.
    ndArray<bool> done;

    //! Lock access to constraint-variables for thread-safety.
    boost::mutex    store_mutex;

    //! Allow access to members for stream-operator
    friend std::ostream & operator<<(std::ostream &os, const OptimizationResult& res);
};

//! Output results as formated human-readable text-stream
std::ostream & operator<<(std::ostream &os, const OptimizationResult& res);

#endif // OPTIMIZATION_RESULT_H
