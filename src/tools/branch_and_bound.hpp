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

#ifndef BRANCH_AND_BOUND_HPP
#define BRANCH_AND_BOUND_HPP

#include <vector>
#include <boost/function.hpp>
#include "thread_pool.hpp"

/*!
 * \brief The BranchAndBound class implements a general branch and bound algorithm.
 *
 * The algorithm is executed on some object of type T by calling a specifiable member-function. The function must return
 * a std::vector of objects. This vector should be empty, if bounding is done for this branch. Otherwise, the returned
 * children are added to a stack of objects which are in turn processed by calling the bound_fcn.
 *
 * The number of threads used in the branch and bound algorithm is user-specifiable. For thread-control the
 * ThreadPool-class is used!
 */
class BranchAndBound
{
public:
    /*!
     * \brief Constructor.
     * \param num_threads num_threads Maximum number of threads to use. If num_threads<2 run() is executed in the
     * calling thread.
     */
    BranchAndBound(size_t num_threads=0)
    {
        if (num_threads > 1) {
           thread_pool = boost::shared_ptr<ThreadPool>(new ThreadPool(num_threads));
        }
    }
    //! Destructor.
    ~BranchAndBound() {}

    /*!
     * \brief Execute branch and bound algorithm.
     * \param root Object on which the bound_fcn should be called.
     * \param bound_fcn Member function which is called in each bound-step. Must return vector of children if further
     * bounding in this direction is necessary.
     *
     * For the single-threaded case this function directly calls bound() with the root-object and member-function as
     * arguments.
     *
     * For the multi-threaded case this function schedules a call to bound() using the root-object and specified member
     * function. and waits until the ThreadPool is done.
     */
    template <class T>
    void run(T& root, std::vector<T> (T::*bound_fcn)()) const
    {
        if (thread_pool.get() == 0) {
            bound(root,bound_fcn);
        } else {
            boost::function<void (void)> f = boost::bind(&BranchAndBound::bound<T>, this, root, bound_fcn);
            thread_pool->schedule(f);
            thread_pool->waitDone();
        }
    }

private:
    /*!
     * \brief Bound-step of branch and bound algorithm.
     * \param node Object on which the bound_fcn should be called.
     * \param bound_fcn Member function which should be called. Must return vector of children if further bounding in
     * this direction is necessary.
     *
     * For the single-threaded case this function recursively calls bound() for all children returned by one call to
     * bound().
     *
     * For the multi-threaded case this function schedules calls to bound() on the ThreadPool for all children returned
     * by one call to bound().
     */
    template <class T>
    void bound(T& node, std::vector<T> (T::*bound_fcn)()) const
    {
        std::vector<T> ret = (node.*bound_fcn)();

        for (size_t i = 0; i<ret.size(); i++) {
            if (thread_pool.get() == 0) {
                 bound(ret[i],bound_fcn);
            } else {
                boost::function<void (void)> f = boost::bind(&BranchAndBound::bound<T>, this, ret[i], bound_fcn);
                thread_pool->schedule(f);
            }
        }
    }

    //! ThreadPool on which the BranchAndBound algorithm is executed for the multi-threaded case.
    boost::shared_ptr<ThreadPool> thread_pool;
};

#endif //BRANCH_AND_BOUND_HPP
