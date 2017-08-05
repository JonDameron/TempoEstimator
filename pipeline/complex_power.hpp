#pragma once

#include <memory>
#include <string>
#include "pipeline/abstract_pipeline_node.hpp"
#include "utility/proc_data.hpp"
#include "utility/util.hpp"

namespace pipeline {

class ComplexPower : public AbstractPipelineNode
{
public:

  ComplexPower ();

  ~ComplexPower ();

private:

  /**
   The # work units should be calculated dynamically based on estimated overall
   computational complexity; for demonstration purposes, though, a reasonable
   constant value is fine.
   */
  static const int kDefaultNWorkUnits = 4;

  std::string Process (std::shared_ptr<const ProcData> input,
                       std::shared_ptr<ProcData>* output) OVERRIDE;
};

}
