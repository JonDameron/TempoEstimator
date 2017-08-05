#pragma once

#include <condition_variable>
#include <memory>
#include <mutex>
#include <string>
#include <vector>
#include <boost/thread/thread.hpp>
#include "pipeline/abstract_pipeline_node.hpp"
#include "utility/proc_data.hpp"
#include "utility/thread_squad.hpp"
#include "utility/util.hpp"

namespace pipeline {

/**
 * Warning: While AbstractPipelineHead is internally thread-safe, it is the
 * responsibility of the external scope to ensure that its public API is
 * accessed in a thread-safe manner.
 */
class AbstractPipelineHead : public AbstractPipelineNode
{
  /** Generally accepted convention is that the 'friend' designation should be
   * used sparingly, but in consideration of the special relationship between
   * the pipeline head and nodes, its use here is tolerated to ease design
   * and avoid sacrificing access appropriate access limitations.
   * Specifically, AbstractPipelineHead needs to invoke
   * AbstractPipelineNode::SetHead(), which really needs to specify access level
   * 'private' to ensure non-pipeline classes cannot accidentally invoke it.
   */
  friend class AbstractPipelineNode;

public:

  explicit AbstractPipelineHead (std::shared_ptr<ThreadSquad> thread_squad);

  virtual ~AbstractPipelineHead ();

  std::string pipeline_error_str () const {
    std::unique_lock<std::mutex> lock (const_cast<std::mutex&>(state_mutex_));
    return pipeline_error_str_;
  }

  std::string Start ();

  /** Stop pipeline execution and block until the main pipeline thread
   * has exited.
   * Calling this method while another thread is waiting for either
   * WaitForProcessingCompletion() or Stop() to return can result in
   * undefined behavior.
   */
  void Stop ();

  /** Block until the main pipeline thread has exited.  Note that the main
   * thread won't exit until all pending thread work units have been
   * fully processed.
   * Calling this method while another thread is waiting for either
   * WaitForProcessingCompletion() or Stop() to return can result in
   * undefined behavior.
   */
  void WaitForProcessingCompletion ();

protected:

  /** Intended to be invoked by CollectHeadData() or indirectly by an
   * implementation of AbstractPipelineNode::Process() via
   * AbstractPipelineNode::InterruptibleSleep(); pauses execution of current
   * thread for given duration in seconds, or until the pipeline enters a
   * shutting-down state, whichever comes first.
   * @returns Reason for return
   */
  ReturnReason InterruptibleSleep (double seconds);

private:

  void MainPipelineThreadFunc ();

  void SubmitThreadWorkToHead (ThreadSquad::ThreadWorkUnit func);

  std::string Process (std::shared_ptr<const ProcData> input_UNUSED,
                       std::shared_ptr<ProcData>* output) OVERRIDE;

  /** Use of interrupt_condvar is optional
   */
  virtual std::string CollectHeadData (std::shared_ptr<ProcData>* data_out) = 0;

  bool is_running_;

  std::shared_ptr<boost::thread> main_thread_;

  std::mutex state_mutex_;

  std::condition_variable state_change_condvar_;

  std::shared_ptr<ThreadSquad> thread_squad_;

  std::string pipeline_error_str_;
};

}
