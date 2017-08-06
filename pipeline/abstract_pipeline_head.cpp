#include <cassert>
#include <functional>
#include <memory>
#include <string>
#include <boost/thread/thread.hpp>
#include "pipeline/abstract_pipeline_head.hpp"
#include "utility/thread_squad.hpp"

using namespace std;
using namespace pipeline;

AbstractPipelineHead :: AbstractPipelineHead (shared_ptr<ThreadSquad> thread_squad)
: pipeline_debug_enabled_(false),
  is_running_(false),
  thread_squad_(thread_squad)
{
}

AbstractPipelineHead :: ~AbstractPipelineHead ()
{
  Stop();
}

string AbstractPipelineHead :: Start ()
{
  if (is_running_) {
    // This should probably be a failed assertion since the caller ought to
    // know better than to invoke Start() twice without an intermediate Stop()
    return "Pipeline has already been started";
  }

  AbstractPipelineNode* node = this;
  while (node)
  {
    node->SetHead(this);
    if (pipeline_debug_enabled_) {
      node->set_node_debug_enabled(true);
    }
    node = node->output_node().get();
  }

  is_running_ = true;
  main_thread_.reset( new boost::thread(
      bind(&AbstractPipelineHead::MainPipelineThreadFunc, this) ) );

  return "";
}

void AbstractPipelineHead :: Stop ()
{
  // Thread synchronization unneeded here since if !is_running_, then the
  // only other thread that accesses is_running_ -- the main thread -- has been
  // halted or is currently halting
  if (!is_running_) {
    return;
  }
  // Sanity check -- if this fails, then the caller is likely using the
  // AbstractPipelineHead public interface in a thread-unsafe manner
  assert(main_thread_);

  // Necessary to ensure thread synchronization since the derived
  // CollectHeadData() can optionally implement an interruptible wait
  { unique_lock<mutex> lock (state_mutex_);
    is_running_ = false;
    state_change_condvar_.notify_all();
  }

  if (main_thread_->joinable()) {
    main_thread_->join();
  }
  main_thread_.reset();

  // Note that it is AbstractPipelineHead::MainPipelineThreadFunc() that
  // pushes data through the pipeline. Thus, now that it's stopped, we can
  // return with the guarantee that the pipeline is completely "stopped".
  // Because the main thread always waits for all work units submitted
  // by each call to AbstractPipelineNode::Process() to be fully evaluated
  // before continuing, we are certain that when MainPipelineThreadFunc()
  // returns, no thread is performing any work for this pipeline.

  AbstractPipelineNode* node = this;
  while (node) {
    node->InvalidateHead();
    node = node->output_node().get();
  }
}

void AbstractPipelineHead :: WaitForProcessingCompletion ()
{
  // MainPipelineThreadFunc exits after the final node of the pipeline has
  // finished processing.  Not checking is_running_ because its state is
  // irrelevant; as long as the main pipeline thread is active, data may still
  // be moving through the pipeline.
  if (main_thread_->joinable()) {
    main_thread_->join();
  }
}

AbstractPipelineNode::ReturnReason AbstractPipelineHead :: InterruptibleSleep (
    double seconds)
{
  unique_lock<mutex> lock (state_mutex_);

  if (!is_running_) {
    return ReturnReason::kShuttingDown;
  }

  cv_status wait_result =
      state_change_condvar_.wait_for(lock, chrono::duration<double>(seconds));

  // If our running state changed without a corresponding signalling of
  // state_change_condvar_, then there's a bug in the AbstractPipelineHead;
  // it should not be possible to change is_running_ to false without this
  // condvar being signalled
  assert(cv_status::no_timeout == wait_result || is_running_);

  return is_running_ ? ReturnReason::kTimeExpired : ReturnReason::kShuttingDown;
}

void AbstractPipelineHead :: MainPipelineThreadFunc ()
{
  const pid_t this_thread_lwpid = GetCurrentThreadId();

  /*
   NOTE: Throughout this function, we're able to safely check is_running_ without
   first acquiring state_mutex_ since the only method that can alter is_running_
   while this main thread is active is Stop(), which will then patiently wait
   on thread::join() before continuing.
   */

  { unique_lock<mutex> lock (state_mutex_);
    pipeline_error_str_.clear();
  }

  while (is_running_)
  {
    shared_ptr<const ProcData> input;
    AbstractPipelineNode* node = this;

    while (node && is_running_)
    {
      string result;
      shared_ptr<ProcData> output;

      // Note that <input> will be null on the initial pass; this is intentional
      result = node->Process(input, &output);

      // Any thread work units submitted by the derived AbstractPipelineNode's
      // Process() implementation must be completed before the output ProcData
      // can travel onward to the next pipeline node
      thread_squad_->WaitForWorkUnits(this_thread_lwpid);

      if (!result.empty()) {
        node->set_process_error_str(result);
        { unique_lock<mutex> lock (state_mutex_);
          pipeline_error_str_ = result;
        }
        return;
      }
      if (!output) {
        // Assume no error occurred, but we've reached the end of the pipeline.
        // Note that, unlike in Stop(), we *don't* need to worry about acquiring
        // state_mutex_ before altering is_running_ here since it's guaranteed
        // that the derived CollectHeadData() -- the only method that could be
        // adversely affected by this action -- isn't invoked, as that method
        // is only called by this thread.
        is_running_ = false;
        break;
      }

      if (OutputPreservationMode::kReference == node->output_preservation_mode()) {
        node->set_most_recent_output(output);
      } else if (OutputPreservationMode::kCopy == node->output_preservation_mode()) {
        node->set_most_recent_output(ProcData::New(*output));
      } else {
        // Don't make any special effort to preserve the latest output, and
        // blow away any previously preserved output
        node->set_most_recent_output(shared_ptr<ProcData>());
      }

      // The next node in the pipeline will receive the current node's output
      // as input
      input = output;

      // There should be no danger of the returned node becoming invalid (instance
      // being freed) while we're using it since every AbstractPipelineNode refers
      // to its output as a shared_ptr; thus, as long as the head node is valid,
      // it's guaranteed that every other node in the pipeline is also valid.
      node = node->output_node().get();
    }
  }
}

void AbstractPipelineHead :: SubmitThreadWorkToHead (
    ThreadSquad::ThreadWorkUnit func)
{
  // This method is private, and will only be invoked by
  // AbstractPipelineNode::Process() implementations of pipeline components.

  // The public ThreadSquad API is internally thread-safe, so we don't need to
  // worry about regulating multi-threaded access to this method.

  // Further, since AbstractPipelineNode::Process() is guaranteed to be invoked
  // only from within MainPipelineThreadFunc(), and this method is invoked only
  // from within the former, the fact that execution has made it to this point
  // implies that the main head func is running, and the pipeline state is
  // therefore valid.
  // In consideration of this, there is no need to acquire a mutex here.
  assert(thread_squad_);

  thread_squad_->PushThreadWork(func);
}

string AbstractPipelineHead :: Process (
    shared_ptr<const ProcData> input_UNUSED, shared_ptr<ProcData>* output)
{
  // Prevent "unused parameter" warning
  (void)input_UNUSED;

  if (!is_running_) {
    // Not setting <output> informs the caller that the input stream
    // has terminated
    return "";
  }

  return CollectHeadData(output);
}
