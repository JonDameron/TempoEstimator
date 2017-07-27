#include <mutex>
#include <boost/thread/thread.hpp>
#include "utility/thread_squad.hpp"

using namespace std;

ThreadSquad :: ThreadSquad (int n_threads)
: is_running_(true)
{
  for (int thread_id = 0; thread_id < n_threads; ++thread_id)
  {
    threads_.emplace_back(bind(&ThreadSquad::MainThreadFunc, this, thread_id));
  }
}

ThreadSquad :: ~ThreadSquad ()
{
  { unique_lock<mutex> lock (squad_state_mutex_);
    work_queue_.clear();
    is_running_ = false;
    squad_state_changed_condvar_.notify_all();
  }

  for (int i = 0; i < threads_.size(); ++i)
  {
    threads_[i].join();
  }
}

void ThreadSquad :: PushThreadWork (ThreadWorkUnit work_unit)
{
  { unique_lock<mutex> lock (squad_state_mutex_);
    work_queue_.push_back(work_unit);
    squad_state_changed_condvar_.notify_one();
  }
}

void ThreadSquad :: MainThreadFunc (int thread_id)
{
  ThreadWorkUnit work_unit;

  while (true)
  {
    { unique_lock<mutex> lock (squad_state_mutex_);
      while (is_running_ && work_queue_.empty()) {
        squad_state_changed_condvar_.wait(lock);
      }
      if (!is_running_) {
        return;
      }
      work_unit = work_queue_.front();
      work_queue_.pop_front();
    }

    string result = work_unit(thread_id);
    if (!result.empty()) {
      // TODO: More sophisticated error handling
      cerr << "Thread " << thread_id << ": Error while executing work unit: "
           << result << "\n";
    }
  }
}
