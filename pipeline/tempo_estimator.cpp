/* Musical Tempo Estimation Library
 * (c) 2017 Jonathan G. Dameron | jon.g.dameron@gmail.com
 */

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include "pipeline/tempo_estimator.hpp"
#include "utility/util.hpp"

using namespace std;
using namespace pipeline;

TempoEstimator :: TempoEstimator (
        const FftProcessorConfig& level_1_fft_cfg,
        const FftProcessorConfig& level_2_fft_cfg,
        std::shared_ptr<ComplexPower> complex_power_pipeline_node,
        double min_acceptable_tempo_samples_per_beat,
        double max_acceptable_tempo_samples_per_beat)
: level_1_fft_cfg_(level_1_fft_cfg),
  level_2_fft_cfg_(level_2_fft_cfg),
  complex_power_pipeline_node_(complex_power_pipeline_node),
  min_acceptable_tempo_samples_per_beat_(min_acceptable_tempo_samples_per_beat),
  max_acceptable_tempo_samples_per_beat_(max_acceptable_tempo_samples_per_beat),
  tempo_offset_samples_(-1),
  average_tempo_samples_per_beat_(0)
{
  // We're not influencing the behavior of complex_power_pipeline_node in any
  // way, we're just mandating that a copy of its output be preserved if the
  // caller hasn't already set this up; note that the caller changing this to
  // OutputPreservationMode::kNone AFTER this constructor returns will likely
  // lead to a failed assertion later.
  if (AbstractPipelineNode::OutputPreservationMode::kNone ==
      complex_power_pipeline_node->output_preservation_mode())
  {
    complex_power_pipeline_node->set_output_preservation_mode(
        AbstractPipelineNode::OutputPreservationMode::kReference);
  }
}

TempoEstimator :: ~TempoEstimator ()
{
}

double TempoEstimator :: EstimateTempoOffsetNSamples (
    double tempo_n_samples_per_beat,
    int level_1_reference_bin_index)
{
  /* Correlate the level-1 power with a short-tick pulse train of frequency
   * corresponding to max_group_peak_bin_index, and note the time offset at
   * which the correlation sum is maximized
   */

  // FIXME: This will crash for *very* short input songs or very long level-1 FFTs

  string res;

  const int level_1_fft_len      = level_1_fft_cfg_.fft_len();
  const int level_1_out_n_bins   = level_1_fft_cfg_.out_len_per_fft();
  const int level_1_n_ffts_total = level_1_fft_cfg_.n_ffts();

  // Number of level-1 FFTs we're actually going to take into consideration
  // when performing the correlation; we're primarily interested in lining up
  // the tempo with an early portion of the song.
  const int level_1_n_ffts_proc = level_1_n_ffts_total / 4;

  const int n_samples_total = level_1_fft_len * level_1_n_ffts_proc;

  const double n_beats_total =
      n_samples_total / tempo_n_samples_per_beat;

  const double level_1_n_ffts_per_beat = level_1_n_ffts_proc / n_beats_total;

  const int corr_n_steps = (int)(level_1_n_ffts_per_beat + 1);

  const int beat_pulse_train_len = level_1_n_ffts_proc - corr_n_steps;

  // Shorter is generally better, as this reduces the worst-case error when
  // attempting to work out the optimal alignment of the beat pulse train later,
  // but making it too short can increase the risk of total detection failure,
  // as the amount of averaging is reduced.
  const int beat_pulse_len = std::max(1, (int)(level_1_n_ffts_per_beat / 8));

  const int n_pulses = (int)(beat_pulse_train_len / level_1_n_ffts_per_beat) - 2;

  const int beat_pulse_offset = 0;

  // Set the 'pulse' portions of our pulse train to 1
  vector<double> beat_pulse_train (beat_pulse_train_len, 0);
  for (int beat_pulse_i = 0; beat_pulse_i < n_pulses; ++beat_pulse_i)
  {
    // We want one pulse per beat
    const int train_offset =
        (int)(beat_pulse_offset + beat_pulse_i * level_1_n_ffts_per_beat);

    std::fill( beat_pulse_train.begin() + train_offset,
               beat_pulse_train.begin() + train_offset + beat_pulse_len,
               double(1) );
  }

  // This condition should have been assured via set_output_preservation_mode()
  // in the TempoEstimator constructor
  assert(complex_power_pipeline_node_->most_recent_output());

  const double* level_1_pwr =
      complex_power_pipeline_node_->most_recent_output()->double_data();

  const double* corr_in_1  = level_1_pwr + level_1_reference_bin_index;
  const int corr_in_1_step = level_1_out_n_bins;
  const int corr_in_1_len  = level_1_n_ffts_proc;
  const double* corr_in_2  = beat_pulse_train.data();
  const int corr_in_2_step = 1;
  const int corr_in_2_len  = beat_pulse_train.size();

  // Sanity check, this should have been calculated properly above
  assert(corr_in_1_len >= corr_in_2_len);

  // The length of the correlation output is input #1 len - input #2 len + 1
  vector<double> corr_output_vec (corr_in_1_len - corr_in_2_len + 1);

  res = Correlate( corr_in_1,
                   corr_in_1_step,
                   corr_in_1_len,
                   corr_in_2,
                   corr_in_2_step,
                   corr_in_2_len,
                   corr_output_vec.data(),
                   corr_output_vec.size() );
  // If Correlate() returned an error, then there's a software bug
  assert(res.empty());

  auto corr_peak_iter = std::max_element(corr_output_vec.begin(),
                                         corr_output_vec.end());

  const int corr_peak_index = corr_peak_iter - corr_output_vec.begin();

  // The correlation peak index properly aligns the beat pulse train over the
  // audio input power vector, but we need to add a fraction of the width of
  // the 'pulse' component of the train to account for the fact that the beat
  // is somewhere within the pulse duration -- probably within the first half
  // of the pulse, as many instruments that would contribute to this detection
  // method have a relatively short 'attack' time and longer 'decay' time.
  const double beat_pulse_alignment_fac = 0.5;
  const double offset_n_samples =
      (corr_peak_index + beat_pulse_len * beat_pulse_alignment_fac)
      * level_1_fft_len;

  return offset_n_samples;
}

string TempoEstimator :: Process (shared_ptr<const ProcData> input,
                                  shared_ptr<ProcData>* output)
{
  // Instead of setting *output, TempoEstimator sets internal variables.
  // The final calculations are accessible via public getters.
  // *output being unset will be interpreted as an indication that this
  // pipeline node is an endpoint.

  // <input> must be the complex-valued output of the second-level FFT.

  /*
   ** Tempo Estimation Algorithm Description **
   We use the phases associated with the highest-amplitude bins of the
   level-2 FFT to approximate the tempo offset (time offset of first beat).
   In the context of this processing module, the term "group" refers to a
   collection of cornerturn level-2 FFTs spanning a full level-1 range of
   output bins (accounting for gaps between level-2 FFTs per
   level_2_fft_cfg.in_dist).
   TODO: Finish this documentation
   */

  const complex<double>* level_2_in = input->complex_double_data();

  // The output of a r2c transform is symmetric, and FFTW omits the mirror
  // image:  http://www.fftw.org/fftw3_doc/Real_002ddata-DFT-Array-Format.html
  const int level_1_out_n_bins = level_1_fft_cfg_.out_len_per_fft();
  const int level_2_out_n_bins = level_2_fft_cfg_.out_len_per_fft();

  vector<double> level_2_in_pwr (input->n_complex_doubles()); // complex => real
  for (int i = 0, i_end = input->n_complex_doubles(); i < i_end; ++i)
  {
    level_2_in_pwr.at(i) = Sqr(level_2_in[i].real()) + Sqr(level_2_in[i].imag());
  }

  assert(0 == level_1_out_n_bins % level_2_fft_cfg_.in_dist());
  const int n_level_2_ffts_per_group =
      level_1_out_n_bins / level_2_fft_cfg_.in_dist();

  const int n_groups = level_1_fft_cfg_.n_ffts() / level_2_fft_cfg_.fft_len();

  vector<double> level_2_in_pwr_group_avg (
      level_2_in_pwr.size() / n_level_2_ffts_per_group );

  // g: group index
  // b: level-2 bin index
  // f: level-2 FFT index
  for (int g = 0; g < n_groups; ++g)
  {
    const int g_start = g * n_level_2_ffts_per_group * level_2_out_n_bins;
    for (int b = 0; b < level_2_out_n_bins; ++b)
    {
      const int b_start = g_start + b;
      double sum = 0;
      for (int f = 0; f < n_level_2_ffts_per_group; ++f)
      {
        sum += level_2_in_pwr.at(b_start + f*level_2_out_n_bins);
      }
      level_2_in_pwr_group_avg.at(g*level_2_out_n_bins + b) =
          sum / n_level_2_ffts_per_group;
    }
  }

  // index of max, val of max
  vector<pair<int, double>> level_2_fft_peaks (level_2_fft_cfg_.n_ffts(),
                                               pair<int, double>(-1, -1));

  // f: level-2 FFT index
  for (int f = 0; f < level_2_fft_cfg_.n_ffts(); ++f)
  {
    auto dc_pwr_bin_iter = level_2_in_pwr.begin() + f*level_2_out_n_bins;
    // Note that we're excluding several bins at/near the output extremities,
    // including the first (DC) bin, as candidates for max element, as
    // these -- especially DC -- can be extremely large and have bin indices
    // too small to contribute to tempo detection.
    auto max_pwr_iter = std::max_element(
        dc_pwr_bin_iter + 2, dc_pwr_bin_iter + level_2_out_n_bins - 1);
    level_2_fft_peaks.at(f).first  = max_pwr_iter - dc_pwr_bin_iter;
    level_2_fft_peaks.at(f).second = *max_pwr_iter;
  }

  vector<pair<double, double>> level_2_weighted_fft_peaks (
      level_2_fft_cfg_.n_ffts(), pair<double, double>(-1, -1) );

  // f: level-2 FFT index
  for (int f = 0; f < level_2_fft_cfg_.n_ffts(); ++f)
  {
    // try to get a more accurate max index by taking a weighted average of
    // the peak and its two adjacent bins

    const int peak_bin = level_2_fft_peaks.at(f).first;
    auto dc_pwr_bin_iter = level_2_in_pwr.begin() + f*level_2_out_n_bins;

    // Use fourth root of the level-2 FFT output power for this operation
    // since sqrt(LevelTwoPower) => LevelTwoMagnitude => LevelOnePower
    // and sqrt(LevelOnePower) => LevelOneMagnitude
    vector<double> bins_mag = {
        pow(*(dc_pwr_bin_iter + peak_bin - 1), 0.25),
        pow(*(dc_pwr_bin_iter + peak_bin    ), 0.25),
        pow(*(dc_pwr_bin_iter + peak_bin + 1), 0.25)
    };
    double bins_mag_sum = std::accumulate<vector<double>::const_iterator, double>(
        bins_mag.begin(), bins_mag.end(), 0 );
    double index_adjust =
        -1 * bins_mag.at(0) / bins_mag_sum + 1 * bins_mag.at(2) / bins_mag_sum;
    level_2_weighted_fft_peaks.at(f).first  = peak_bin + index_adjust;
    // pow(4) because we took the fourth root earlier
    level_2_weighted_fft_peaks.at(f).second = pow(bins_mag_sum/3, 4);
  }

  // index of max, val of max
  vector<pair<int, double>> level_2_fft_group_peaks (n_groups,
                                                     pair<int, double>(-1, -1));

  // g: group index
  for (int g = 0; g < n_groups; ++g)
  {
    auto dc_pwr_bin_iter =
        level_2_in_pwr_group_avg.begin() + g*level_2_out_n_bins;

    // Note that we're excluding several bins at/near the output extremities,
    // including the first (DC) bin, as candidates for max element, as
    // these -- especially DC -- can be extremely large and have bin indices
    // too small to contribute to tempo detection.
    auto max_pwr_iter = std::max_element(
        dc_pwr_bin_iter + 2, dc_pwr_bin_iter + level_2_out_n_bins - 1);

    level_2_fft_group_peaks.at(g).first  = max_pwr_iter - dc_pwr_bin_iter;
    level_2_fft_group_peaks.at(g).second = *max_pwr_iter;
  }

  vector<pair<double, double>> level_2_weighted_fft_group_peaks (
      n_groups, pair<double, double>(-1, -1) );

  // g: group index
  for (int g = 0; g < n_groups; ++g)
  {
    // Try to get a more accurate max index by taking a weighted average of
    // the peak and its two adjacent bins

    const int group_peak_bin = level_2_fft_group_peaks.at(g).first;
    auto dc_pwr_bin_iter =
        level_2_in_pwr_group_avg.begin() + g*level_2_out_n_bins;

    // Use fourth root of the level-2 FFT output power for this operation
    // since sqrt(LevelTwoPower) => LevelTwoMagnitude => LevelOnePower
    // and sqrt(LevelOnePower) => LevelOneMagnitude
    vector<double> peak_bins_mag = {
        pow(*(dc_pwr_bin_iter + group_peak_bin - 1), 0.25),
        pow(*(dc_pwr_bin_iter + group_peak_bin    ), 0.25),
        pow(*(dc_pwr_bin_iter + group_peak_bin + 1), 0.25)
    };
    // Note the importance of explicitly casting the accumulator's initial value
    // to double; this ensures that the numeric type of the sum is 'double'
    const double peak_bins_mag_sum =
        std::accumulate(peak_bins_mag.begin(), peak_bins_mag.end(), double(0));
    const double peak_index_adjust =
        (peak_bins_mag.at(2) - peak_bins_mag.at(0)) / peak_bins_mag_sum;

    level_2_weighted_fft_group_peaks.at(g).first =
        group_peak_bin + peak_index_adjust;
    // pow(4) because we took the fourth root earlier
    level_2_weighted_fft_group_peaks.at(g).second = pow(peak_bins_mag_sum/3, 4);
  }

  auto weighted_peak_compare_func =
      [](const pair<double, double>& v1, const pair<double, double>& v2)->bool {
        return v1.second < v2.second;
      };

  auto max_group_peak_bin_iter =
      std::max_element( level_2_weighted_fft_group_peaks.begin(),
                        level_2_weighted_fft_group_peaks.end(),
                        weighted_peak_compare_func );

  const double max_group_peak_bin_index = max_group_peak_bin_iter->first;

#if 0
  // Dump level-2 power data to stderr as CSV for postproc/testing in Matlab.
  // If this code is enabled, BE SURE to redirect stderr to a file
  // (or elsewhere), otherwise the terminal will be mercilessly inundated with
  // an onslaught of comma-separated values.
  // Example: tempo_estimator -i HueyLewisBackInTime.wav 2> you_bet_marty.csv
  std::for_each( level_2_in_pwr.begin(), level_2_in_pwr.end(),
      [level_2_out_n_bins, &level_2_in_pwr](const double& pwr) {
        const int i = &pwr - level_2_in_pwr.data();
        cerr << pwr << (level_2_out_n_bins-1 == i % level_2_out_n_bins ? "\n" : ",");
      } );
#endif

  // TODO: Eventually only print the following diagnostic info if a 'verbose'
  // option is enabled

  cout << "\n";

  cout << "Level-1 FFT configuration:\n" << level_1_fft_cfg_.ToString() << "\n\n";
  cout << "Level-2 FFT configuration:\n" << level_2_fft_cfg_.ToString() << "\n\n";

#if 1
  for (int f = 0; f < level_2_fft_cfg_.n_ffts(); ++f)
  {
    cout << "level_2_weighted_fft_peaks[" << f << "] = "
         << level_2_weighted_fft_peaks.at(f).first
         << " , " << level_2_weighted_fft_peaks.at(f).second << "\n";
  }
#endif

  cout << "\n";

  for (int g = 0; g < n_groups; ++g)
  {
    cout << "level_2_weighted_fft_group_peaks[" << g << "] = "
         << level_2_weighted_fft_group_peaks.at(g).first
         << " , " << level_2_weighted_fft_group_peaks.at(g).second
         << "; level_2_fft_group_peaks[" << g << "] = "
         << level_2_fft_group_peaks.at(g).first
         << " , " << level_2_fft_group_peaks.at(g).second
         << "\n";
  }

  // As implied by the names, some of these values are counts of time-domain
  // audio samples.

  const int n_time_samples_per_level_1_fft = level_1_fft_cfg_.fft_len();

  const int max_beat_rate_bin_index = level_2_out_n_bins - 1;

  // Number of samples spanning the minimum beat interval (i.e., time between
  // beats at the maximum beat rate).
  // This is essentially our way of compensating for an odd-length level-2 FFT,
  // the max evaluated frequency of which has a period slightly greater than
  // two samples.
  const double n_time_samples_per_max_beat_rate_interval =
      level_2_fft_cfg_.fft_len() / (double)max_beat_rate_bin_index
      * n_time_samples_per_level_1_fft;

  double tempo_estimate_samples_per_beat =     // Nyquist => max measurable tempo
      n_time_samples_per_max_beat_rate_interval * (double)max_beat_rate_bin_index
      / max_group_peak_bin_index;

  if (tempo_estimate_samples_per_beat < min_acceptable_tempo_samples_per_beat_)
  {
    // Scale it up by the smallest power of two that would place it above the
    // 'acceptable' minimum. Using a power of two so we won't get a polyrhythmic
    // effect for 4/4 time signature music.
    const int linear_scale = 1 + (int)(min_acceptable_tempo_samples_per_beat_
                                       / tempo_estimate_samples_per_beat);
    const int exp_scale = 1 << (int)ceil(log2(linear_scale));

    tempo_estimate_samples_per_beat *= exp_scale;
  }
  else if (tempo_estimate_samples_per_beat > max_acceptable_tempo_samples_per_beat_)
  {
    // Scale it down by the smallest integer that would place it below the
    // 'acceptable' maximum
    const int scale = 1 + (int)(tempo_estimate_samples_per_beat
                                / max_acceptable_tempo_samples_per_beat_);

    tempo_estimate_samples_per_beat /= scale;
  }

  const int group_start_bin_for_tempo_offset_ref = 0; // Use first group

  const int level_2_check_bin_for_tempo_offset_ref =
      (int)(0.5 + max_group_peak_bin_index);

  const int level_2_fft_start_index_for_tempo_offet_ref = 0;

  double level_2_max_pwr_for_tempo_offset_ref = -1;
  int level_2_max_pwr_fft_index_for_tempo_offset_ref = -1;
  // f: level-2 FFT index
  for (int f = level_2_fft_start_index_for_tempo_offet_ref;
       f < n_level_2_ffts_per_group; ++f)
  {
    const int check_bin =
        group_start_bin_for_tempo_offset_ref
        + f * level_2_out_n_bins
        + level_2_check_bin_for_tempo_offset_ref;

    const double pwr = level_2_in_pwr.at(check_bin);
    if (pwr > level_2_max_pwr_for_tempo_offset_ref) {
      level_2_max_pwr_for_tempo_offset_ref = pwr;
      level_2_max_pwr_fft_index_for_tempo_offset_ref = f;
    }
  }

  // What we're trying to determine here is which level-1 frequency bin would
  // the best choice when correlating the level-1 output powers with a beat
  // pulse train of pulse frequency matching our estimated tempo
  const int level_1_reference_bin_index_for_tempo_offset_estimation =
        level_2_max_pwr_fft_index_for_tempo_offset_ref
        * level_2_fft_cfg_.in_dist();

  tempo_offset_samples_ = EstimateTempoOffsetNSamples(
      tempo_estimate_samples_per_beat,
      level_1_reference_bin_index_for_tempo_offset_estimation );

  average_tempo_samples_per_beat_ = tempo_estimate_samples_per_beat;

  if (!isfinite(tempo_offset_samples_) || tempo_offset_samples_ < 0) {
    return "Failed to calculate tempo offset, result is '"
        + to_string(tempo_offset_samples_) + "'";
  }
  if (!isfinite(average_tempo_samples_per_beat_)
        || average_tempo_samples_per_beat_ <= 0) {
    return "Failed to calculate average tempo, result is '"
        + to_string(average_tempo_samples_per_beat_) + "'";
  }

  /*
  TODO: Remove this; mapping phase to its corresponding time offset yields
        inaccurate results for most songs.
  const int in_index_for_tempo_offset =
      level_2_out_n_bins * level_2_fft_cfg_.n_ffts()/2
      + level_2_fft_group_peaks.at(0).first;

  const double positive_tempo_bin_phase =
      fmod(2*M_PI + std::arg(level_2_in[in_index_for_tempo_offset]), 2*M_PI);

  assert(in_index_for_tempo_offset < input->n_complex_doubles());
  tempo_offset_samples_ = 1/(2*M_PI) * average_tempo_samples_per_beat_ *
      positive_tempo_bin_phase;
  */

  return "";
}
