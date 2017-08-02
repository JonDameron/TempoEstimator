#include <iostream>
#include <mutex>
#include <boost/thread/thread.hpp>
#include "utility/thread_squad.hpp"

using namespace std;

ThreadSquad :: ThreadSquad (int n_threads)
: is_running_(true)
{
  for (int thread_id = 0; thread_id < n_threads; ++thread_id)
  {
    threads_.emplace_back(bind(&ThreadSquad::ThreadFunc, this, thread_id));
  }
}

ThreadSquad :: ~ThreadSquad ()
{
  { unique_lock<mutex> lock (squad_state_mutex_);
    work_queue_.clear();
    is_running_ = false;
    squad_state_changed_condvar_.notify_all();

    for (WorkCountInfo& work_count_info : work_count_by_submitter_) {
      unique_lock<mutex> work_count_info_lock (work_count_info.mux);
      work_count_info.count_is_zero_condvar.notify_all();
    }
    work_count_by_submitter_.clear();
  }

  for (int i = 0; i < threads_.size(); ++i)
  {
    threads_[i].join();
  }
}

void ThreadSquad :: PushThreadWork (ThreadWorkUnit work_unit)
{
  unique_lock<mutex> squad_state_lock (squad_state_mutex_);

  work_queue_.emplace_back(GetCurrentThreadId(), work_unit);
  squad_state_changed_condvar_.notify_one();

  // GetCurrentThreadId() returns a unique ID (lightweight PID) for each thread
  const pid_t this_thread_lwpid = GetCurrentThreadId();
  pair<pid_t, WorkCountInfo> new_map_item (this_thread_lwpid, WorkCountInfo(0));

  // If the current thread already has its own mapping in
  // work_count_by_submitter_, then the returned insert() iterator
  // (insert_result.first) will point to it instead of a new entry created
  // for new_map_item
  auto insert_result = work_count_by_submitter_.insert(new_map_item);

  // insert_result.first is an iterator to the item of interest, which is
  // a pair<KeyType, ValueType>.
  // insert_result.second is a bool indicating whether an insertion was necessary.

  // It's safe to lock a WorkCountInfo mutex while also possessing the squad
  // state mutex
  unique_lock<mutex> work_count_info_lock (insert_result.first->second.mux);
  ++insert_result.first->second.count;
}

void ThreadSquad :: WaitForWorkUnits (pid_t thread_id_of_submitter)
{
  WorkCountInfo* info = nullptr;

  { unique_lock<mutex> squad_state_lock (squad_state_mutex_);
    auto map_iter = work_count_by_submitter_.find(thread_id_of_submitter);
    if (work_count_by_submitter_.end() == map_iter) {
      return; // no work is pending for given thread ID (LWPID)
    }
    info = &map_iter->second;
    // Essential to acquire this mutex before relinquishing squad_state_mutex_
    info->mux.lock();
  }

  // Note that we already have ownership of info->mux
  { unique_lock<mutex> work_count_info_lock (info->mux, std::adopt_lock_t());
    if (0 == info->count) {
      return; // no work is pending for given thread ID (LWPID)
    }
    info->count_is_zero_condvar.wait(work_count_info_lock);
  }
}

void ThreadSquad :: ThreadFunc (int thread_id)
{
  while (true)
  {
    QueuedWorkUnit queued_work;

    { unique_lock<mutex> lock (squad_state_mutex_);
      while (is_running_ && work_queue_.empty()) {
        squad_state_changed_condvar_.wait(lock);
      }
      if (!is_running_) {
        return;
      }
      queued_work = work_queue_.front();
      work_queue_.pop_front();
      // NOTE: Not decrementing the corresponding WorkCountInfo::count until the
      // work unit has actually been executed
    }

    // Execute the work unit function
    string result = queued_work.work_unit(thread_id);

    // Regardless of execution result, decrement pending work unit counter for
    // the thread that submitted the work
    { unique_lock<mutex> lock (squad_state_mutex_);
      auto work_count_info_iter = work_count_by_submitter_.find(queued_work.submitter_pid);
      // There are conceivable cases, e.g. ThreadSquad is undergoing destruction,
      // in which this map find() operation can fail
      if (work_count_by_submitter_.end() != work_count_info_iter)
      {
        // It's safe to lock a WorkCountInfo mutex while also possessing the
        // squad state mutex
        unique_lock<mutex> work_count_info_lock (work_count_info_iter->second.mux);
        --work_count_info_iter->second.count;
        // Sanity check
        assert(work_count_info_iter->second.count >= 0);
        if (0 == work_count_info_iter->second.count) {
          work_count_info_iter->second.count_is_zero_condvar.notify_all();
          // TODO: Eventually remove mapped WorkCountInfo instances when count
          // reaches zero to save a small amount of memory
        }
      }
    }

    if (!result.empty()) {
      // TODO: More sophisticated error handling
      cerr << "Thread " << thread_id << ": Error while executing work unit: "
           << result << "\n";
    }
  }
}
