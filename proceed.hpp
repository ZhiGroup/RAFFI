/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 * This file is responsible for synchronization.
 *
 */

#ifndef PROCEED_HPP
#define PROCEED_HPP

#include <mutex>
#include <condition_variable>
#include <vector>

class Proceed {
public:
    std::mutex mutex;

    Proceed(int num_threads);

    int increment_num_blocked_and_get();

    void update_num_blocked();

    bool can_thread_proceed(int thread);

    void allow_thread_proceed(int thread);

    void disallow_thread_proceed(int thread);

    void allow_all_threads_proceed();

    bool has_all_blocked();

    void mark_thread_finished(int thread);

    bool has_all_finished();

    int get_num_finished();

    void wait_for_master(std::unique_lock<std::mutex> &lock);

    void wait_for_workers(std::unique_lock<std::mutex> &lock);

    void signal_master();

    void signal_workers();

private:
    int num_threads;
    int num_blocked = 0;
    int num_finished = 0;
    std::vector<bool> can_proceed;
    std::vector<bool> have_finished;
    std::condition_variable master_wait_on_this;
    std::condition_variable worker_wait_on_this;
};

#endif


