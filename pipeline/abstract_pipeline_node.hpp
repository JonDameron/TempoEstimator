#pragma once

#include <memory>
#include <string>
#include "utility/proc_data.hpp"
#include "pipeline/abstract_pipeline_head.hpp"

namespace pipeline {

class AbstractPipelineNode
{
public:

  AbstractPipelineNode ();

  virtual ~AbstractPipelineNode ();

  std::string ConnectOutput (std::shared_ptr<AbstractPipelineNode> output_node);

protected:

  /** Returns an empty string on success or a string describing the error.
   */
  typedef std::function<std::string()> ThreadWorkFunc;

  std::shared_ptr<AbstractPipelineNode> output_node () {
    return output_node_;
  }

  void SubmitThreadWork (ThreadWorkFunc func);

private:

  virtual std::string Process (std::shared_ptr<const ProcData> input,
                               std::shared_ptr<ProcData>* output) = 0;

  std::shared_ptr<AbstractPipelineHead> head_;
  std::shared_ptr<AbstractPipelineNode> output_node_;
};

}
