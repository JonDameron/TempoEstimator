#pragma once

#include <memory>
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

  shared_ptr<ThreadSquad> thread_squad_;

  std::string Process (std::shared_ptr<const ProcData> input_UNUSED,
                       std::shared_ptr<ProcData>* output) override;

  virtual std::string CollectHeadData (std::shared_ptr<ProcData>* data_out) = 0;
};

}
