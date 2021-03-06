#pragma once

#include <condition_variable>
#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
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
   * thread_index is the unique zero-based index of the thread -- NOT its LWPID.
   * See also GetCurrentThreadId(), which returns a thread's LWPID.
   */
  typedef std::function<std::string (int thread_index)> ThreadWorkUnit;

  /** Construct a new ThreadSquad of given max number of concurrent threads.
   */
  explicit ThreadSquad (int n_threads);

  /** Destructor; clears the work queue and stops all running threads.
   */
  ~ThreadSquad ();

  /** Push a new work unit to the back of the queue (lowest priority).
   */
  void PushThreadWork (ThreadWorkUnit work_unit);

  void WaitForWorkUnits (pid_t thread_id_of_submitter);

private:

  void ThreadFunc (int thread_id);

  bool is_running_;

  struct QueuedWorkUnit
  {
    QueuedWorkUnit ()
    : submitter_pid(0)
    {}
    QueuedWorkUnit (pid_t submitter_pid_in, ThreadWorkUnit work_unit_in)
    : submitter_pid(submitter_pid_in), work_unit(work_unit_in)
    {}
    pid_t submitter_pid;
    ThreadWorkUnit work_unit;
  };
  std::deque<QueuedWorkUnit> work_queue_;

  std::vector<boost::thread> threads_;

  std::mutex squad_state_mutex_;

  std::condition_variable squad_state_changed_condvar_;

  /** Number of pending work units organized by the PID of the submitting thread.
   * Note that we avoid over-utilizing the general squad state mutex/condvar by
   * assigning a unique mutex/condvar to each work count.  Since the original
   * motivation for creating this mapping is to enable a thread to wait for
   * any given submitter's work units to fully evaluate, basic thread sync
   * operations are necessary.
   */
  struct WorkCountInfo
  {
    WorkCountInfo (int initial_count)
    : count(initial_count),
      mux(new std::mutex()),
      count_is_zero_condvar(new std::condition_variable())
    {}
    int count;
    // NOTE: It would be lovely to skip the shared_ptr wrapping here and instead
    // simply use unordered_map::emplace with piecewise_construct_t to insert
    // a new WorkCountInfo in-place into work_count_by_submitter_, but
    // unfortunately gcc < 4.8.0 doesn't support the emplace() method.
    // Bottom line, for gcc < 4.8.0, we need a copy-constructible WorkCountInfo.
    std::shared_ptr<std::mutex> mux;
    std::shared_ptr<std::condition_variable> count_is_zero_condvar;
  };
  std::unordered_map<pid_t, WorkCountInfo> work_count_by_submitter_;
};
