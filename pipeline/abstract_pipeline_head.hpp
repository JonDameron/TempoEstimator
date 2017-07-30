#pragma once

#include <condition_variable>
#include <memory>
#include <mutex>
#include <string>
#include <vector>
#include <boost/thread/thread.hpp>
#include "utility/proc_data.hpp"
#include "utility/thread_squad.hpp"
#include "pipeline/abstract_pipeline_node.hpp"

namespace pipeline {

class AbstractPipelineHead : public AbstractPipelineNode
{
public:

  explicit AbstractPipelineHead (std::shared_ptr<ThreadSquad> thread_squad);

  virtual ~AbstractPipelineHead ();

  std::string Start ();

  void Stop ();

private:

  bool is_running_;

  std::shared_ptr<boost::thread> main_thread_;

  std::mutex state_mutex_;

  std::condition_variable state_change_condvar_;

  std::shared_ptr<ThreadSquad> thread_squad_;

  std::string MainPipelineThreadFunc ();

  std::string Process (std::shared_ptr<const ProcData> input_UNUSED,
                       std::shared_ptr<ProcData>* output) override;

  /** Use of interrupt_condvar is optional
   */
  virtual std::string CollectHeadData (std::condition_variable& interrupt_condvar,
                                       std::shared_ptr<ProcData>* data_out) = 0;
};

}
