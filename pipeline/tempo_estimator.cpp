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
  const int beat_pulse_len = std::max(1, (int)(level_1_n_ffts_per_beat * 0.3)); // * 0.125));

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
  vector<double> corr_out (corr_in_1_len - corr_in_2_len + 1);

  res = Correlate( corr_in_1,
                   corr_in_1_step,
                   corr_in_1_len,
                   corr_in_2,
                   corr_in_2_step,
                   corr_in_2_len,
                   corr_out.data(),
                   corr_out.size() );
  // If Correlate() returned an error, then there's a software bug
  assert(res.empty());

  auto corr_peak_iter = std::max_element(corr_out.begin(),
                                         corr_out.end());

  const int corr_peak_index = corr_peak_iter - corr_out.begin();

  // The correlation peak index properly aligns the beat pulse train over the
  // audio input power vector, but we need to add a fraction of the width of
  // the 'pulse' component of the train to account for the fact that the beat
  // is somewhere within the pulse duration -- probably within the first half
  // of the pulse, as many instruments that would contribute to this detection
  // method have a relatively short 'attack' time and longer 'decay' time.
  const double beat_pulse_alignment_fac = 0.1; //0.25; //0.5;
  const double offset_n_samples =
      (corr_peak_index + beat_pulse_len * beat_pulse_alignment_fac)
      * level_1_fft_len;

  return offset_n_samples;
}

double TempoEstimator :: CalcLevelTwoTempoBinUsingPopularityMethod (
    const vector<pair<double, double>>& level_2_weighted_fft_peaks)
{
  /* The input level_2_weighted_fft_peaks matrix dimensions are
   * nGroups X nFfts. Each value is a pair with first = *weighted* index of the
   * level-2 FFT output bin whose corresponding power is the max for that FFT,
   * and second = power of bin corresponding to index <first>.
   * We're going to set up a new array of length level2FftOutLen that transforms
   * the input as follows:
   * init newarr [ 0..level2FftOutLen-1 ] = 0
   * (this sum deliberately ignores power in favor of # appearances)
   * newarr [ floor(0.5 + in[grp, ffti].first) ] += 1
   *
   * Then, we correlate the new array with a 2-point boxcar vector to identify
   * the 'most popular' tempo detection bin.
   *
   * Next, sort the input by weighted index and locate the start of the cluster
   * of interest corresponding to the result of the boxcar correlation, move to
   * the approximate center of the cluster, and move in both directions until the
   * difference in adjacent weighted indices is sufficiently large (indicating
   * departure from said cluster).
   *
   * Finally, perform a weighted average of the 'popular' cluster by multiplying
   * fourth root of power by weighted index, summing, and dividing the result
   * by the sum of the fourth roots of the powers. The final weighted average
   * is the floating-point bin index to be returned to the application following
   * some basic conversions.
   *
   * Fourth root explanation:
   * sqrt(LevelTwoPower) => LevelTwoMagnitude => LevelOnePower
   * and sqrt(LevelOnePower) => LevelOneMagnitude
   */

  const int n_groups = level_1_fft_cfg_.n_ffts() / level_2_fft_cfg_.fft_len();

  const int n_level_2_ffts_per_group =
      level_1_fft_cfg_.out_len_per_fft() / level_2_fft_cfg_.in_dist();

  // Verify that the caller is invoking this method correctly
  assert(level_2_weighted_fft_peaks.size() ==
         n_groups * n_level_2_ffts_per_group);

  string res;

  // Adding 1 to length to account for theoretically possible weighting
  // correction placing the new index just past the end.
  // The element type is 'double' because the vec will be passed to the
  // correlation function.
  // Note that we're ignoring the bins near the FFT output extremeties.

  vector<double> popularity_vec (level_2_fft_cfg_.out_len_per_fft() + 1, 0);
  const int popularity_vec_min_bin =
      std::max(6, (int)(0.002 * popularity_vec.size()));
  const int popularity_vec_max_bin =
      (int)(0.997 * popularity_vec.size());
  for (const pair<double, double>& weighted_peak : level_2_weighted_fft_peaks)
  {
    const int peak_bin_index = (int)(0.5 + weighted_peak.first);
    if (peak_bin_index >= popularity_vec_min_bin
          && peak_bin_index <= popularity_vec_max_bin) {
      ++ popularity_vec.at(peak_bin_index);
    }
  }

  vector<double> boxcar (2, 0.5);

  const double* corr_in_1  = popularity_vec.data();
  const int corr_in_1_step = 1;
  const int corr_in_1_len  = popularity_vec.size();
  const double* corr_in_2  = boxcar.data();
  const int corr_in_2_step = 1;
  const int corr_in_2_len  = boxcar.size();

  vector<double> corr_out (corr_in_1_len - corr_in_2_len + 1);

  res = Correlate( corr_in_1,
                   corr_in_1_step,
                   corr_in_1_len,
                   corr_in_2,
                   corr_in_2_step,
                   corr_in_2_len,
                   corr_out.data(),
                   corr_out.size() );
  // If Correlate() returned an error, then there's a software bug
  assert(res.empty());

  auto popularity_begin_iter =
      std::max_element(corr_out.begin(), corr_out.end());
  const int popularity_begin_bin_index = popularity_begin_iter - corr_out.begin();

  vector<pair<double, double>> sorted_peaks = level_2_weighted_fft_peaks;
  std::sort( sorted_peaks.begin(), sorted_peaks.end(),
      [](const pair<double, double>& v1, const pair<double, double>& v2)->bool {
          return v1.first < v2.first;
        } );

  const double popularity_sum_begin_index = popularity_begin_bin_index;

  auto peak_begin_iter = std::find_if(
      sorted_peaks.begin(), sorted_peaks.end(),
      [popularity_sum_begin_index](const pair<double, double>& val)->bool {
          return (int)(0.5 + val.first) >= popularity_sum_begin_index;
        } );

  // Move backward until the step distance exceeds a threshold (i.e., until
  // we seem to have exited the cluster of interest)

  // FIXME: This was originally 0.5 but has to be 1 until I generalize the
  // correlation input vector to have a granularity matching this value;
  // see the granularity parameter for the 'narrow-cluster' calculations below.
  // For now 1 should be fine since the narrow-cluster processing that follows
  // will isolate densely packed peak groups and focus exclusively on them.
  const double cluster_contiguity_threshold = 1.0;

  // Condition "> begin()" which skips the first element is intentional
  auto wide_cluster_begin_iter = peak_begin_iter;
  for (; wide_cluster_begin_iter > sorted_peaks.begin(); --wide_cluster_begin_iter)
  {
    if (wide_cluster_begin_iter->first - (wide_cluster_begin_iter-1)->first
            > cluster_contiguity_threshold) {
      break;
    }
  }

  auto wide_cluster_end_iter = sorted_peaks.end();
  // Sanity check -- this condition would break the following loop, and should
  // never happen
  assert(wide_cluster_begin_iter != sorted_peaks.end());
  for (auto iter = wide_cluster_begin_iter; true; ++iter)
  {
    // Putting the loop break condition in the body just because otherwise the
    // code is really ugly
    if (iter == sorted_peaks.end()-1
          || (iter+1)->first - iter->first > cluster_contiguity_threshold) {
      wide_cluster_end_iter = iter+1;
      break;
    }
  }

  const int wide_cluster_len = wide_cluster_end_iter - wide_cluster_begin_iter;
  const double wide_cluster_peak_range =
      (wide_cluster_end_iter-1)->first - wide_cluster_begin_iter->first;

  // TODO: Really need to organize and document this code. This function has
  // grown out of control; originally it was much shorter.

  // Narrow our peak bin range selection even further by performing another
  // even shorter boxcar correlation on our initial, wide range of bins

  // Make a copy of the 'wide' cluster scaled up by a relatively large factor;
  // the goal is to perform a correlation on a 'zoomed-in' version of our wide
  // cluster, to identify the most dense sub-cluster therein

  const double cluster_zoom_granularity = 1e4;

  // We don't need to include the peak weight in cluster_zoom_vec, just the
  // peak index
  vector<double> cluster_zoom_vec (wide_cluster_len, 0);
  std::transform(
      wide_cluster_begin_iter,
      wide_cluster_end_iter,
      cluster_zoom_vec.begin(),
      [cluster_zoom_granularity](const pair<double, double>& val) {
            return val.first * cluster_zoom_granularity;
          }
      );

  // Adding +1 to account for rounding that could possibly make the size of this
  // vector one value too short
  const int cluster_zoom_pop_vec_len =
      (int)(0.5 + wide_cluster_peak_range * cluster_zoom_granularity) + 1;
  vector<double> cluster_zoom_popularity_vec (cluster_zoom_pop_vec_len, 0);

  for (double cluster_zoom_val : cluster_zoom_vec)
  {
    // Remember that cluster_zoom_vec is sorted since its source vector,
    // sorted_peaks, is sorted
    const int wide_pop_vec_index =
        (int)(0.5 + cluster_zoom_val - cluster_zoom_vec.front());
    ++ cluster_zoom_popularity_vec.at(wide_pop_vec_index);
  }

  const int cluster_zoom_boxcar_len =
      std::min( cluster_zoom_pop_vec_len,
                std::max(2, (int)(0.1 * cluster_zoom_pop_vec_len)) );

  // Note that the sum of the weights of our correlation vector *must* equal the
  // length of the vector (non-normalized) since downstream code assumes that
  // the correlation output is an unaveraged sum.
  // TODO: Try using a triangular shape + pedestal instead of a boxcar so that
  //       outliers from the narrow cluster will have less influence, and the
  //       correlation peak will tend more toward the center of the max overlap.
  //       When using the triangular vector, a special variant of the standard
  //       discrete correlation should be employed that initially positions the
  //       second vector (our triangle) at index -secondVecLen/2 and ends the
  //       correlation at index firstVecLen, of course only performing
  //       multiplications from valid indices 0 -> firstVecLen-1, and
  //       compensating for the secondVec positions at which the vectors are
  //       only partially overlapped by scaling up the result by a factor of
  //       sum(secondVec) / sum(overlappingSecondVec).

  vector<double> cluster_zoom_boxcar (cluster_zoom_boxcar_len, 1);

  const double* cluster_zoom_corr_in_1  = cluster_zoom_popularity_vec.data();
  const int cluster_zoom_corr_in_1_step = 1;
  const int cluster_zoom_corr_in_1_len  = cluster_zoom_popularity_vec.size();
  const double* cluster_zoom_corr_in_2  = cluster_zoom_boxcar.data();
  const int cluster_zoom_corr_in_2_step = 1;
  const int cluster_zoom_corr_in_2_len  = cluster_zoom_boxcar.size();

  vector<double> cluster_zoom_corr_out (
      cluster_zoom_corr_in_1_len - cluster_zoom_corr_in_2_len + 1 );

  res = Correlate( cluster_zoom_corr_in_1,
                   cluster_zoom_corr_in_1_step,
                   cluster_zoom_corr_in_1_len,
                   cluster_zoom_corr_in_2,
                   cluster_zoom_corr_in_2_step,
                   cluster_zoom_corr_in_2_len,
                   cluster_zoom_corr_out.data(),
                   cluster_zoom_corr_out.size() );
  // If Correlate() returned an error, then there's a software bug
  assert(res.empty());

  auto cluster_zoom_corr_peak_begin_iter =
      std::max_element(cluster_zoom_corr_out.begin(), cluster_zoom_corr_out.end());

  const int cluster_zoom_corr_peak_begin_index =
      cluster_zoom_corr_peak_begin_iter - cluster_zoom_corr_out.begin();

  const int cluster_zoom_corr_peak_len = *cluster_zoom_corr_peak_begin_iter;

  const int min_acceptable_cluster_zoom_len = 3;

  vector<pair<double, double>>::iterator weighted_sum_input_begin_iter;
  vector<pair<double, double>>::iterator weighted_sum_input_end_iter;

  if (cluster_zoom_corr_peak_len >= min_acceptable_cluster_zoom_len)
  {
    const double cluster_zoom_popularity_sum_begin_index =
        wide_cluster_begin_iter->first
        + cluster_zoom_corr_peak_begin_index / cluster_zoom_granularity;

    auto narrow_cluster_begin_iter = std::find_if(
        wide_cluster_begin_iter, wide_cluster_end_iter,
        [cluster_zoom_popularity_sum_begin_index](const pair<double, double>& val)->bool {
            return val.first >= cluster_zoom_popularity_sum_begin_index;
          } );

    assert(narrow_cluster_begin_iter != wide_cluster_end_iter);

    const int wide_cluster_begin_index_zoomed_offset =
        (cluster_zoom_corr_peak_begin_iter - cluster_zoom_corr_out.begin())
        + cluster_zoom_boxcar.size();

    const double cluster_zoom_popularity_sum_end_index =
        wide_cluster_begin_iter->first
        + wide_cluster_begin_index_zoomed_offset / cluster_zoom_granularity;

    auto narrow_cluster_end_iter = std::find_if(
        narrow_cluster_begin_iter, wide_cluster_end_iter,
        [cluster_zoom_popularity_sum_end_index](const pair<double, double>& val)->bool {
            return val.first > cluster_zoom_popularity_sum_end_index;
          } );

    assert(narrow_cluster_end_iter > narrow_cluster_begin_iter);

    weighted_sum_input_begin_iter = narrow_cluster_begin_iter;
    weighted_sum_input_end_iter   = narrow_cluster_end_iter;
  }
  else
  {
    weighted_sum_input_begin_iter = wide_cluster_begin_iter;
    weighted_sum_input_end_iter   = wide_cluster_end_iter;
  }

  // Perform weighted sum (iter->first is the peak bin index, second is
  // the weight)
  double weighted_sum = 0;
  double sum_of_weights = 0;
  for (auto iter = weighted_sum_input_begin_iter;
       iter < weighted_sum_input_end_iter; ++iter)
  {
    // To save time, overwrite weights of our (function-local) vector with their
    // fourth roots (see explanation at the top of this function)
    iter->second = pow(iter->second, 0.25);
    weighted_sum += iter->first * iter->second;
    sum_of_weights += iter->second;
  }

  // Calculate average by dividing by the sum of the weights
  const double tempo_bin_relative_to_level_2_fft = weighted_sum / sum_of_weights;

  return tempo_bin_relative_to_level_2_fft;
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
      level_2_in_pwr.size() / n_level_2_ffts_per_group,
      1.0 );

  const int level_2_fft_start = (int)(n_level_2_ffts_per_group * 0.2);
  const int level_2_fft_stop  = (int)(n_level_2_ffts_per_group * 0.8);

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
      for (int f = level_2_fft_start; f < level_2_fft_stop; ++f)
      {
        sum += level_2_in_pwr.at(b_start + f*level_2_out_n_bins);
      }
      level_2_in_pwr_group_avg.at(g*level_2_out_n_bins + b) =
          sum / (level_2_fft_stop - level_2_fft_start);
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
    // TODO: Provide formal documentation for the final values of the many
    //       constants once I'm done experimenting
    auto max_pwr_iter = std::max_element(
        dc_pwr_bin_iter + (int)(level_2_out_n_bins * 0.05),
        dc_pwr_bin_iter + (int)(level_2_out_n_bins * 0.40) ); //0.95) );
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
    const double bins_mag_sum =
        std::accumulate<vector<double>::const_iterator, double>(
            bins_mag.begin(), bins_mag.end(), 0 );
    // For better clarity, this index_adjust calculation could also be written as
    //  -1 * bins_mag.at(0) / bins_mag_sum
    // + 0 * bins_mag.at(1) / bins_mag_sum
    // + 1 * bins_mag.at(2) / bins_mag_sum
    const double index_adjust = (bins_mag.at(2) - bins_mag.at(0)) / bins_mag_sum;
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
    // For better clarity, this index_adjust calculation could also be written as
    //  -1 * peak_bins_mag.at(0) / peak_bins_mag_sum
    // + 0 * peak_bins_mag.at(1) / peak_bins_mag_sum
    // + 1 * peak_bins_mag.at(2) / peak_bins_mag_sum
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

  if (node_debug_enabled())
  {
    cout << "\n";
    cout << "Level-1 FFT configuration:\n" << level_1_fft_cfg_.ToString() << "\n\n";
    cout << "Level-2 FFT configuration:\n" << level_2_fft_cfg_.ToString() << "\n\n";

    for (int f = 0; f < level_2_fft_cfg_.n_ffts(); ++f)
    {
      cout << "level_2_weighted_fft_peaks[" << f << "] = "
           << level_2_weighted_fft_peaks.at(f).first
           << " , " << level_2_weighted_fft_peaks.at(f).second << "\n";
    }

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

    cout << "\n";
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

  // TODO: Finalize this when I'm done experimenting
  //const double level_2_bin_index_for_tempo_estimate = max_group_peak_bin_index;
  const double level_2_bin_index_for_tempo_estimate =
      CalcLevelTwoTempoBinUsingPopularityMethod(level_2_weighted_fft_peaks);

  double tempo_estimate_samples_per_beat =     // Nyquist => max measurable tempo
      n_time_samples_per_max_beat_rate_interval * (double)max_beat_rate_bin_index
      / level_2_bin_index_for_tempo_estimate;

  if (tempo_estimate_samples_per_beat < min_acceptable_tempo_samples_per_beat_)
  {
    // Scale it up by the smallest power of two that would place it above the
    // 'acceptable' minimum. Using a power of two so we won't get a polyrhythmic
    // effect for 4/4 time signature songs.
    const int linear_scale = 1 + (int)(min_acceptable_tempo_samples_per_beat_
                                       / tempo_estimate_samples_per_beat);
    const int exp_scale = 1 << (int)ceil(log2(linear_scale));

    tempo_estimate_samples_per_beat *= exp_scale;
  }
  else if (tempo_estimate_samples_per_beat > max_acceptable_tempo_samples_per_beat_)
  {
    // Scale it down by the smallest power of two that would place it below the
    // 'acceptable' maximum. Using a power of two so we won't get a polyrhythmic
    // effect for 4/4 time signature songs.
    const int linear_scale = 1 + (int)(tempo_estimate_samples_per_beat
                                       / max_acceptable_tempo_samples_per_beat_);
    const int exp_scale = 1 << (int)ceil(log2(linear_scale));

    tempo_estimate_samples_per_beat /= exp_scale;
  }

  // Use second group if it exists
  const int group_begin_bin_for_tempo_offset_ref =
      std::min(1, n_groups-1) * n_level_2_ffts_per_group * level_2_out_n_bins;

  const int level_2_check_bin_for_tempo_offset_ref =
      (int)(0.5 + max_group_peak_bin_index);

  const int level_2_fft_begin_index_for_tempo_offet_ref = 0;

  double level_2_max_pwr_for_tempo_offset_ref = -1;
  int level_2_max_pwr_fft_index_for_tempo_offset_ref = -1;
  // f: level-2 FFT index
  for (int f = level_2_fft_begin_index_for_tempo_offet_ref;
       f < n_level_2_ffts_per_group; ++f)
  {
    const int check_bin =
        group_begin_bin_for_tempo_offset_ref
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
