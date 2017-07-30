#include <functional>
#include <memory>
#include <string>
#include <boost/thread/thread.hpp>
#include "pipeline/abstract_pipeline_head.hpp"
#include "utility/thread_squad.hpp"

using namespace std;
using namespace pipeline;

AbstractPipelineHead :: AbstractPipelineHead (shared_ptr<ThreadSquad> thread_squad)
: is_running_(false),
  thread_squad_(thread_squad)
{
}

AbstractPipelineHead :: ~AbstractPipelineHead ()
{
}

string AbstractPipelineHead :: Start ()
{
  if (is_running_) {
    return "Main pipeline thread is already running";
  }

  is_running_ = true;
  main_thread_.reset( new boost::thread(
      bind(&AbstractPipelineHead::MainPipelineThreadFunc, this) ) );
  return "";
}

void AbstractPipelineHead :: Stop ()
{
  // Necessary to ensure thread synchronization since the derived
  // CollectHeadData() can optionally implement an interruptible wait
  // via state_change_condvar_
  { unique_lock<mutex> lock (state_mutex_);
    is_running_ = false;
    state_change_condvar_.notify_one();
  }

  main_thread_->join();
}

string AbstractPipelineHead :: MainPipelineThreadFunc ()
{
  shared_ptr<const ProcData> input;
  AbstractPipelineNode* node = this;

  while (node)
  {
    string result;
    shared_ptr<ProcData> output;

    { unique_lock<mutex> lock (state_mutex_);
      if (!is_running_) {
        break;
      }
      // Note that <input> will initially be null; this is intentional
      result = node->Process(input, &output);
    }
    if (!result.empty()) {
      stringstream strm (result);
      return StackError(&strm, "Pipeline head failed to collect new data");
    }
    if (!output) {
      // Assume no error occurred, but we've reached the end of the input stream
      is_running_ = false;
      break;
    }

    // The next node in the pipeline will receive this node's output as input
    input = output;

    // There should be no danger of the returned node becoming invalid (instance
    // being freed) while we're using it since every AbstractPipelineNode refers
    // to its output as a shared_ptr; thus, as long as the head node is valid,
    // it's guaranteed that every other node in the pipeline is also valid.
    node = node->output_node().get();
  }

  return "";
}

string AbstractPipelineHead :: Process (
    shared_ptr<const ProcData> input_UNUSED, shared_ptr<ProcData>* output)
{
  // Prevent "unused parameter" warning
  (void)input_UNUSED;

  // Note that the caller, which is guaranteed to be MainPipelineThreadFunc(),
  // is responsible for acquiring state_mutex_

  return CollectHeadData(state_change_condvar_, output);
}
