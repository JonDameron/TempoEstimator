#pragma once

#include <memory>
#include <string>
#include <fftw3.h>
#include "pipeline/abstract_pipeline_node.hpp"
#include "pipeline/fft_processor_config.hpp"
#include "utility/proc_data.hpp"
#include "utility/util.hpp"

namespace pipeline {

class FftProcessor : public AbstractPipelineNode
{
public:

  explicit FftProcessor (const FftProcessorConfig& cfg);

  ~FftProcessor ();

private:

  /**
   The # work units should be calculated dynamically based on estimated overall
   computational complexity; for demonstration purposes, though, a reasonable
   constant value is fine.
   */
  static const int kDefaultNWorkUnits = 4;

  std::string InitPlans (int n_work_units);

  std::string Process (std::shared_ptr<const ProcData> input,
                       std::shared_ptr<ProcData>* output) OVERRIDE;

  FftProcessorConfig cfg_;
  std::vector<fftw_plan> fftw_plans_;
  std::vector<int> work_unit_fft_offsets_;
  std::vector<int> work_unit_n_ffts_;
};

}
