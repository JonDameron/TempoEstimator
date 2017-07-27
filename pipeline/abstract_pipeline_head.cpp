#include <memory>
#include <string>
#include "pipeline/abstract_pipeline_head.hpp"

using namespace std;
using namespace pipeline;

AbstractPipelineHead :: AbstractPipelineHead (shared_ptr<ThreadSquad> thread_squad)
: thread_squad_(thread_squad)
{

}

AbstractPipelineHead :: ~AbstractPipelineHead ()
{

}

string AbstractPipelineHead :: Start ()
{
  return "";
}

void AbstractPipelineHead :: Stop ()
{

}

string AbstractPipelineHead :: Process (
    shared_ptr<const ProcData> input_UNUSED, shared_ptr<ProcData>* output)
{
  // Prevent "unused parameter" warning
  (void)input_UNUSED;

  return CollectHeadData(output);
}
