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

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <boost/shared_ptr.hpp>

#include "function.h"
#include "optimization_result.h"
#include "parameters.h"
#include "vec.hpp"

/*!
 * \brief The class Optimization is either used to optimize a cost function or
 * to solve a set of inequality constraints.
 *
 * It must be noted that the optimize() function is currently not implemented so
 * this class is currently only used for constraint solving. In addition to
 * "normal" inequality constraints the special case of stability analyis using
 * a value set can be treated by passing a Function constructed from a ValueSet.
 *
 * Constraints are solved by passing a Function to constrain(). This uses the
 * BranchAndBound class to execute the algorithm in doConstrain().
 */
class Optimization
{
public:
    //! Construct an instance and initialize \a num_threads to given value.
    Optimization(size_t num_threads);
    //! Destructor.
    virtual ~Optimization() {}

    /*!
     * \brief Determine minimum of Function \a cost using doOptimize(). Not
     * usefull as doOptimize() is currently not implemented!
     */
    void optimize(const Function &cost);
    /*!
     * \brief Determine paramters for which the constraints given in Function
     * constraint are satisfied and stores the results in Parameters::res.
     * \param constraint Function specifying constraints. Can be a constructed
     * from an Equation, i.e., a number of equations or a ValueSet.
     *
     * Executes a BranchAndBound algorithm which calls doConstrain() in each
     * step. The number of threads used when executing the BranchAndBound
     * algorithm is specified in the constructor Optimization().
     *
     * Blocks until the BranchAndBound algorithm is done.
     */
    void constrain(const Function &constraint);

    //! Retrieve OptimizationResult Optimization::res
    OptimizationResult* result() {return res.get();}

private:
    //! Optimization step in BranchAndBound algorithm. CURRENTLY NOT IMPLEMENTED!
    std::vector<Optimization> doOptimize();
    /*!
     * \brief Evaluate constraints, store results, and split parameter set if
     * necessary.
     *
     * Implements the algorithm presented in \cite Schauss2017 when executed
     * within a BranchAndBound algorithm. First performs the "bound"-step.
     * Depending on the result the Function may be split into two functions on
     * smaller parameter sets which are then returned to the BranchAndBound
     * algorithm. This subdivision is performed in splitGradB(), splitGradTm(),
     * or splitPar(), depending on the evaluation result of the Function.
     */
    std::vector<Optimization> doConstrain();

    /*!
     * \brief Subdivide in direction of the maximum estimated partial derivative
     * of the Bernstein polynomial.
     *
     * The subdivision direction is determined by calling Function::dpoly_dparam().
     * Then, splitAt() is called with this direction.
     */
    std::vector<Optimization> splitGradB() const;
    /*!
     * \brief Subdivide in direction of the maximum partial derivative of the
     * interval remainder.
     *
     * The subdivision direction is determined by calling Function::drest_dparam().
     * Note that some handling of special cases is then done, e.g., we do not
     * split in \em sigma-direction for boundary mapping.
     * Finally, splitAt() is called with this direction.
     *
     * \todo In this case this is rather inefficient as the TaylorModels were
     * already evaluated for the different directions. Also, splitAt() actually
     * only performs a Bernstein subdivision while copying the interval.
     * remainder. So a reevaluation will probably be triggered in the next step.
     */
    std::vector<Optimization> splitGradTm() const;
    /*!
     * \brief Subdivide in direction of the maximum relative parameter size.
     *
     * The subdivision direction is determined by calling Parameters::maxRelFractionDir().
     * Then, splitAt() is called with this direction.
     */
    std::vector<Optimization> splitPar() const;

    /*!
     * \brief Subdivide Function in \a direction at \a mid_point by calling
     * Function::splitAt(). Returns vector of Optimization objects.
     *
     * These Optimization objects are initialized to the same value as this
     * Optimization object which means all of the shared pointers, e.g., for
     * statistics are copied. Then, the functions in Optimization::fcn are
     * replaced by the subdivided functionsm Optimization::depth is incremented
     * by one, and a distinct Optimization::id is set for both objects while
     * incrementing Optimization::id_count.
     */
    std::vector<Optimization> splitAt(size_t direction, double mid_point = 0.5) const;

    //! Print progress to std::cout indicating \a state as part of the output.
    void printProgress(const std::string &state) const;

    //! Shared pointer to Function.
    boost::shared_ptr<Function>   fcn;

    //! Shared pointer to OptimizationResult.
    boost::shared_ptr<OptimizationResult> res;

    //! Number of treads to use for optimization or constraint solving set in constructor.
    size_t num_threads;

    //! Number of instantiated Optimization objects.
    boost::shared_ptr<size_t> id_count;
    //! Number of calls to splitB()
    boost::shared_ptr<size_t> split_B_count;
    //! Number of calls to splitTm()
    boost::shared_ptr<size_t> split_Tm_count;
    //! Number of calls to splitPar()
    boost::shared_ptr<size_t> split_par_count;
    /*!
     * \brief Progress of constraint solving or optimization. Should be one
     * when done.
     *
     * \todo Not correct for boundary mapping in value sets as we are actually
     * mapping two value sets, therefore a value of nearly 200% is reached.
     * \todo Does not work currently for stability check of disjoint regions.
     * Sould however by possible to evaluate how much of the right-half plane
     * has been checked.
     */
    boost::shared_ptr<double> fraction_done;

    //! ID of object corresponding to id_count at instantiation.
    size_t id;
    //! Depth of this object in BranchAndBound-tree (starting at zero)
    size_t depth;

    //! Stream operator for debug-output
    friend std::ostream & operator<<(std::ostream &os, const Optimization& opt);
};

//! Stream operator for debug-output. Prints different counts.
std::ostream & operator<<(std::ostream &os, const Optimization& opt);

#endif // OPTIMIZATION_H
