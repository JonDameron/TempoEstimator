#pragma once

#include <condition_variable>
#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <boost/thread/thread.hpp>
#include "utility/util.hpp"

/** Encapsulates a group of threads to execute arbitrary work units.
 */
class ThreadSquad
{
public:

  /** A thread work unit is implemented as a std::function, thereby enabling
   * closure functionality.
   */
  typedef std::function<std::string (int thread_id)> ThreadWorkUnit;

  explicit ThreadSquad (int n_threads);

  ~ThreadSquad ();

  void PushThreadWork (ThreadWorkUnit work_unit);

private:

  void MainThreadFunc (int thread_id);

  bool is_running_;

  std::vector<boost::thread> threads_;

  std::mutex squad_state_mutex_;

  std::condition_variable squad_state_changed_condvar_;

  std::deque<ThreadWorkUnit> work_queue_;
};
