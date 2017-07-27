#include "pipeline/abstract_pipeline_node.hpp"

using namespace std;
using namespace pipeline;

AbstractPipelineNode :: AbstractPipelineNode ()
{
}

AbstractPipelineNode :: ~AbstractPipelineNode ()
{
}

string AbstractPipelineNode :: ConnectOutput (
    shared_ptr<AbstractPipelineNode> output_node)
{
  if (output_node_) {
    return "Pipeline node's output is already connected to another node";
  }

  output_node_ = output_node;

  return "";
}
