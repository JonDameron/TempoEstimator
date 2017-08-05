#include "pipeline/abstract_pipeline_head.hpp"
#include "pipeline/abstract_pipeline_node.hpp"

using namespace std;
using namespace pipeline;

AbstractPipelineNode :: AbstractPipelineNode ()
: head_(nullptr),
  output_preservation_mode_(OutputPreservationMode::kNone)
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

void AbstractPipelineNode :: SubmitThreadWork (ThreadSquad::ThreadWorkUnit work)
{
  AbstractPipelineHead* typed_head = dynamic_cast<AbstractPipelineHead*>(head_);

  // The pipeline head was initially received as an AbstractPipelineHead, so
  // the dynamic cast should always succeed
  assert(typed_head);

  typed_head->SubmitThreadWorkToHead(work);
}

AbstractPipelineNode::ReturnReason AbstractPipelineNode :: InterruptibleSleep (
    double seconds)
{
  AbstractPipelineHead* typed_head = dynamic_cast<AbstractPipelineHead*>(head_);

  // The pipeline head was initially received as an AbstractPipelineHead, so
  // the dynamic cast should always succeed
  assert(typed_head);

  return typed_head->InterruptibleSleep(seconds);
}

void AbstractPipelineNode :: SetHead (AbstractPipelineHead* head)
{
  // This method is for internal use only, so the validity of <head> should
  // be guaranteed
  assert(head);

  head_ = head;
}

void AbstractPipelineNode :: InvalidateHead ()
{
  head_ = nullptr;
}
