/* Musical Tempo Estimation Library
 * (c) 2017 Jonathan G. Dameron | jon.g.dameron@gmail.com
 */

#pragma once

#include <memory>
#include <string>
#include "pipeline/abstract_pipeline_node.hpp"
#include "pipeline/complex_power.hpp"
#include "pipeline/fft_processor.hpp"
#include "utility/pcm_mono_audio_data.hpp"
#include "utility/proc_data.hpp"
#include "utility/util.hpp"

namespace pipeline {

class TempoEstimator : public AbstractPipelineNode
{
public:

  TempoEstimator (const FftProcessorConfig& level_1_fft_cfg,
                  const FftProcessorConfig& level_2_fft_cfg,
                  std::shared_ptr<ComplexPower> complex_power_pipeline_node,
                  double min_acceptable_tempo_samples_per_beat,
                  double max_acceptable_tempo_samples_per_beat);

  ~TempoEstimator ();

  /** Approximate offset, in floating-point # audio samples, of the first 'beat'.
   * The value returned must be divided by the sample rate (waveform
   * samples-per-second) to convert to time.
   * Updated by Process(), invalid until Process() has been invoked at least once.
   */
  double tempo_offset_samples () const {
    return tempo_offset_samples_;
  }

  /** Approximate average musical tempo. Updated by Process().
   * Updated by Process(), invalid until Process() has been invoked at least once.
   */
  double average_tempo_samples_per_beat () const {
    return average_tempo_samples_per_beat_;
  }

private:

  /**
   The # work units should be calculated dynamically based on estimated overall
   computational complexity; for demonstration purposes, though, a reasonable
   constant value is fine.
   */
  static const int kDefaultNWorkUnits = 4;

  double EstimateTempoOffsetNSamples (double tempo_n_samples_per_beat,
                                      int level_1_reference_bin_index);

  double CalcLevelTwoTempoBinUsingPopularityMethod (
      const std::vector<std::pair<double, double>>& level_2_weighted_fft_peaks);

  std::string Process (std::shared_ptr<const ProcData> input,
                       std::shared_ptr<ProcData>* output) OVERRIDE;

  FftProcessorConfig level_1_fft_cfg_;

  FftProcessorConfig level_2_fft_cfg_;

  /** The complex power of the level-1 output spectra is needed for time offset
   * estimation of the first beat (the tempo offset).
   * This is safe since the pipeline only processes a single block of raw input
   * at a time.
   */
  std::shared_ptr<ComplexPower> complex_power_pipeline_node_;

  double min_acceptable_tempo_samples_per_beat_;

  double max_acceptable_tempo_samples_per_beat_;

  double tempo_offset_samples_;

  double average_tempo_samples_per_beat_;
};

}
