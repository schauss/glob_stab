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

#ifndef THREADPOOL_HPP
#define THREADPOOL_HPP

#include <boost/shared_ptr.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/function.hpp>

#include <deque>
#include <vector>

/*!
 * \brief The ThreadPool class implements a simple thread pool based on Boost Thread.
 *
 * It allows to specify the number of threads which should be executed concurrently and offers two different interfaces
 * over which a function call can be scheduled to be executed on the thread pool.
 */
class ThreadPool
{
public:
    /*!
     * \brief Constructor.
     * \param num_threads Number of threads to instantiate and therefore maximum number of jobs which are executed in
     * parallel. Defaults to number of processor cores.
     * \param fifo First-In-First-Out (true) or Last-In-First-Out (false).
     */
    ThreadPool(size_t num_threads = boost::thread::hardware_concurrency(), bool fifo = true)
        : num_threads(num_threads), num_threads_working(0), stop_working(false), fifo(fifo)
    {
		if (num_threads == 0)
			throw("ThreadPool::ThreadPool num_threads==0");

        for (size_t i=0; i<num_threads; i++)
            threads.push_back(boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&ThreadPool::work, this))));
    }

    /*!
     * \brief Destructor.
     *
     * Waits until all currently running and scheduled jobs are done and then destroys the thread pool. Note that
     * calling the destructor does also not prevent new jobs from being scheduled on the thread pool and these will also
     * be executed.
     */
    ~ThreadPool()
    {
        if (!done())
                std::cout << "~ThreadPool: not done -> will wait!" << std::endl;

        stop_working = true;

        for (size_t i=0; i<num_threads; i++)
            threads[i]->join();
    }

    /*!
     * \brief Schedule job to be executed.
     * \param object boost::function to execute.
     */
    void schedule(boost::function<void (void)> object)
    {
        boost::mutex::scoped_lock l(objects_mutex);
        objects.push_back(object);
    }

    /*!
     * \brief Schedule job to be executed.
     * \param object Object of which a member function should be executed.
     * \param mem_fn Member function of object which should be executed.
     */
    template <class T>
    void schedule(const T &object, void (T::*mem_fn)())
    {
        boost::function<void (void)> f = boost::bind(mem_fn, object);
        boost::mutex::scoped_lock l(objects_mutex);
        objects.push_back(f);
    }

    //! Wait until all jobs are done
    void waitDone()
    {
        while (!done())
            usleep(10000);
    }

    //! Are all jobs done?
    bool done()
    {
        boost::mutex::scoped_lock l(objects_mutex);
        return (objects.empty() && (num_threads_working == 0) );
    }

    //! Return number of instantiated threads.
    size_t numThreads() const {return num_threads;}
    //! Return number of threads which are currently working
    size_t numThreadsWorking() const {return num_threads_working;}

private:
    /*!
     * \brief Main function which is executed by each thread.
     *
     * Checks whether there are any jobs scheduled. If there are then pop one job from objects and run it. If there
     * aren't then sleep and recheck periodically.
     */
    void work()
    {
        boost::function<void (void)> my_obj;
        while (!stop_working || !done()) {
            objects_mutex.lock();
            // check if there is an object on stack
            if (!objects.empty()) {

                // take top object from stack
                if (fifo) {
                    my_obj = objects.front();
                    objects.pop_front();
                } else {
                    my_obj = objects.back();
                    objects.pop_back();
                }

                // update number of working threads
                num_threads_working++;
                objects_mutex.unlock();

                // do work
                try {
                    my_obj();
                } catch (std::exception &err) {
                    std::cerr << "std::exception in threadpool: " << err.what() << std::endl;
					throw;
                } catch (const char* err) {
                    std::cerr << "Error in threadpool: " << err << std::endl;
					throw;
                }

                // update number of working threads
                objects_mutex.lock();
                num_threads_working--;
                objects_mutex.unlock();
            } else {
                objects_mutex.unlock();
                usleep(1000);
            }
        }
    }

    //! Double-ended queue of boost::function objects which should be executed (i.e., scheduled jobs).
    std::deque<boost::function<void (void)> > objects;
    boost::mutex objects_mutex; //!< Mutex to lock manipulation of objects queue.

    std::vector<boost::shared_ptr<boost::thread> > threads; //!< Threads in the thread pool.

    size_t num_threads; //!< Number of instantiated threads.
    size_t num_threads_working;  //!< Number of threads which are currently working.

    bool stop_working; //! If true, stop workers as soon as there is no more work left.
    bool fifo; //! Execute jobs FIFO (true) or LIFO (false)
};

#endif // THREADPOOL_HPP
