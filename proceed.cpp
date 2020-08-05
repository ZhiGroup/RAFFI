/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 * This file is responsible for synchronization.
 *
 */

#include <algorithm>

#include "proceed.hpp"

/**
 * Constructor of Proceed.
 *
 * @param num_threads number of worker threads
 */
Proceed::Proceed(int num_threads) :
    num_threads(num_threads),
    can_proceed(num_threads, false),
    have_finished(num_threads, false) {}


/**
 * Increment num_blocked by one and return it.
 */
int
Proceed::increment_num_blocked_and_get()
{
    ++num_blocked;
    return num_blocked;
}


/**
 * Set num_blocked to the number of threads that have not finished.
 */
void
Proceed::update_num_blocked()
{
    num_blocked = std::count(have_finished.begin(), have_finished.end(), true);
}


/**
 * @param thread thread ID.
 *
 * @return whether the thread can resume.
 */
bool
Proceed::can_thread_proceed(int thread)
{
    return can_proceed[thread];
}


/**
 * Allow the input thread to resume.
 *
 * @param thread thread ID.
 */
void
Proceed::allow_thread_proceed(int thread)
{
    can_proceed[thread] = true;
}


/**
 * Disallow the input thread to resume.
 *
 * @param thread thread ID.
 */
void
Proceed::disallow_thread_proceed(int thread)
{
    can_proceed[thread] = false;
}


/**
 * Allow all threads to resume.
 */
void
Proceed::allow_all_threads_proceed()
{
    std::fill(can_proceed.begin(), can_proceed.end(), true);
}


/**
 * @return whether all threads are blocked. Finished threads are considered blocked.
 */
bool
Proceed::has_all_blocked() 
{
    return num_blocked == num_threads;
}


/**
 * Finished threads are considered blocked.
 *
 * @param thread thread ID.
 */
void
Proceed::mark_thread_finished(int thread)
{
    ++num_finished;
    have_finished[thread] = true;
}


/**
 * @return whether all threads have finished.
 */
bool
Proceed::has_all_finished()
{
    return num_finished == num_threads;
}


/**
 * @return number of threads that have finished.
 */
int
Proceed::get_num_finished()
{
    return num_finished;
}


/**
 * Block calling worker thread util resumed by the master thread.
 *
 * @param lock a unique lock of this Proceed.
 */
void
Proceed::wait_for_master(std::unique_lock<std::mutex> &lock) {
    worker_wait_on_this.wait(lock);
}


/**
 * Block master util resumed by the last blocked worker thread.
 *
 * @param lock a unique lock of this Proceed.
 */
void
Proceed::wait_for_workers(std::unique_lock<std::mutex> &lock) {
    master_wait_on_this.wait(lock);
}


/**
 * Resume master thread.
 */
void
Proceed::signal_master() {
    master_wait_on_this.notify_one();
}

/**
 * Resume all worker threads.
 */
void
Proceed::signal_workers() {
    worker_wait_on_this.notify_all();
}
