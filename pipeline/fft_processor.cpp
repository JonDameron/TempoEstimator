#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <fftw3.h>
#include "utility/util.hpp"
#include "pipeline/fft_processor.hpp"
#include "pipeline/fft_processor_config.hpp"

using namespace std;
using namespace pipeline;

FftProcessor :: FftProcessor (const FftProcessorConfig& cfg)
: cfg_(cfg)
{
  InitPlans(kDefaultNWorkUnits);
}

FftProcessor :: ~FftProcessor ()
{
  for (fftw_plan& plan : fftw_plans_) {
    fftw_destroy_plan(plan);
  }
}

string FftProcessor :: InitPlans (int n_work_units)
{
  if (!fftw_plans_.empty()) {
    return "Already initialized";
  }

  const int  rank         = 1;
  const int  fft_len      = cfg_.fft_len();
  const int  n_ffts_total = cfg_.n_ffts();
  const int* in_n_embed   = nullptr;
  const int  in_stride    = cfg_.in_stride();
  const int  in_dist      = cfg_.in_dist();
  const int* out_n_embed  = nullptr;
  const int  out_stride   = cfg_.out_stride();
  const int  out_dist     = cfg_.out_dist();
  const int  sign         = cfg_.sign();
  const int  flags        = cfg_.fftw_planner_flags();
  // applicable only for real-to-real (r2r) transforms
  const fftw_r2r_kind real_to_real_kind = cfg_.real_to_real_kind();

  const vector<int> fft_lengths_for_planner (1, fft_len);
  const vector<fftw_r2r_kind>
      real_to_real_kinds_for_planner (1, real_to_real_kind);

  if (FftProcessorConfig::kInvalidFftwPlannerTimeLimit !=
        cfg_.fftw_planner_time_limit_sec())
  {
    // Warning: fftw_set_timelimit() sets a global parameter. Thus, multiple
    // FftProcessor instances with valid time limit settings
    // (i.e., != kInvalidFftwPlannerTimeLimit) will override one another. The
    // application must be careful to avoid this.
    static int fftw_set_timelimit_count = 0;
    // Print at most one warning.
    if (2 == ++fftw_set_timelimit_count) {
      cerr << "Warning: Multiple valid settings for fftw_planner_time_limit_sec"
           << " override one another\n";
    }
    fftw_set_timelimit(cfg_.fftw_planner_time_limit_sec());
  }

  // These temporary buffers are solely for plan construction; when actually
  // executing the FFT, we'll provide different buffers.
  void* temp_in_array =
      fftw_malloc(cfg_.in_value_size() * cfg_.total_input_span_n_values());
  void* temp_out_array =
      fftw_malloc(cfg_.out_value_size() * cfg_.total_output_span_n_values());

  int fft_offset_this_unit = 0;

  for (int unit_i = 0; unit_i < n_work_units; ++unit_i)
  {
    // If n_ffts_total doesn't evenly divide n_work_units, we'll need to account
    // for the remainder by distributing it across that many work units.
    const int n_ffts_this_unit = n_ffts_total / n_work_units
        + (unit_i < n_ffts_total % n_work_units ? 1 : 0);

    switch (cfg_.type())
    {
    case FftProcessorConfig::FftType::kComplexToComplex:
      fftw_plans_.push_back(
          fftw_plan_many_dft(
              rank,
              fft_lengths_for_planner.data(),
              n_ffts_this_unit,
              static_cast<fftw_complex*>(temp_in_array),
              in_n_embed,
              in_stride,
              in_dist,
              static_cast<fftw_complex*>(temp_out_array),
              out_n_embed,
              out_stride,
              out_dist,
              sign,
              flags )
          );
      break;
    case FftProcessorConfig::FftType::kRealToComplex:
      fftw_plans_.push_back(
          fftw_plan_many_dft_r2c(
              rank,
              fft_lengths_for_planner.data(),
              n_ffts_this_unit,
              static_cast<double*>(temp_in_array),
              in_n_embed,
              in_stride,
              in_dist,
              static_cast<fftw_complex*>(temp_out_array),
              out_n_embed,
              out_stride,
              out_dist,
              flags )
          );
      break;
    case FftProcessorConfig::FftType::kComplexToReal:
      fftw_plans_.push_back(
          fftw_plan_many_dft_c2r(
              rank,
              fft_lengths_for_planner.data(),
              n_ffts_this_unit,
              static_cast<fftw_complex*>(temp_in_array),
              in_n_embed,
              in_stride,
              in_dist,
              static_cast<double*>(temp_out_array),
              out_n_embed,
              out_stride,
              out_dist,
              flags )
          );
      break;
    case FftProcessorConfig::FftType::kRealToReal:
      fftw_plans_.push_back(
          fftw_plan_many_r2r(
              rank,
              fft_lengths_for_planner.data(),
              n_ffts_this_unit,
              static_cast<double*>(temp_in_array),
              in_n_embed,
              in_stride,
              in_dist,
              static_cast<double*>(temp_out_array),
              out_n_embed,
              out_stride,
              out_dist,
              real_to_real_kinds_for_planner.data(),
              flags )
          );
      break;
    }

    work_unit_fft_offsets_.push_back(fft_offset_this_unit);
    work_unit_n_ffts_.push_back(n_ffts_this_unit);

    fft_offset_this_unit += n_ffts_this_unit;
  }

  return "";
}

string FftProcessor :: Process (shared_ptr<const ProcData> input,
                                shared_ptr<ProcData>* output)
{
  stringstream res;

  const int n_work_units = fftw_plans_.size();

  const size_t input_value_size  = cfg_.in_value_size();
  const size_t output_value_size = cfg_.out_value_size();

  const int total_input_span_n_vals  = cfg_.total_input_span_n_values();
  const int total_output_span_n_vals = cfg_.total_output_span_n_values();

  const size_t total_input_span_bytes =
      total_input_span_n_vals * input_value_size;

  if (input->size_bytes() < total_input_span_bytes) {
    res << "Input buffer too short; minimum size = " << total_input_span_bytes
        << ", actual size = " << input->size_bytes();
    return res.str();
  }

  const size_t total_output_span_bytes =
      total_output_span_n_vals * output_value_size;

  *output = ProcData::New(total_output_span_bytes);
  // ProcData::New() should always succeed
  assert(*output);

  for (int unit_i = 0; unit_i < n_work_units; ++unit_i)
  {
    // fftw_plan is an opaque pointer type, so we can just copy it
    fftw_plan plan       = fftw_plans_.at(unit_i);
    const int fft_offset = work_unit_fft_offsets_.at(unit_i);
    const int n_ffts     = work_unit_n_ffts_.at(unit_i);

    // Have to cast away const-ness of the input data since the FFTW API
    // requires a mutable input buffer. Note that despite this, FFTW will not
    // modify the input because the transform plan was set up as out-of-place
    // (different input/output buffers).
    // Further, we're using untyped buffers here since the value of cfg_.type()
    // will determine the I/O types.

    void* input_start = const_cast<uint8_t*>(input->byte_data())
        + input_value_size * fft_offset * cfg_.in_dist();

    void* output_start = (*output)->byte_data()
        + output_value_size * fft_offset * cfg_.out_dist();

    ThreadSquad::ThreadWorkUnit thread_work;

    // IMPORTANT: Because we're using the "new-array" FFTW execution functions
    // here, it is essential that we are careful about our selections for
    // the input and output buffers. As they point to locations within
    // existing larger buffers, we must ensure that (1) the buffer size from
    // <start> to the end of the allocation is sufficient, and (2) the ALIGNMENT
    // is 16-byte, i.e. the given offset is on a 16-byte boundary, for SIMD
    // optimization requirements.

    switch (cfg_.type())
    {
    case FftProcessorConfig::FftType::kComplexToComplex:
      thread_work = [plan, input_start, output_start]()->string {
          fftw_complex* typed_input  = static_cast<fftw_complex*>(input_start);
          fftw_complex* typed_output = static_cast<fftw_complex*>(output_start);
          fftw_execute_dft(plan, typed_input, typed_output);
        };
      break;
    case FftProcessorConfig::FftType::kRealToComplex:
      thread_work = [plan, input_start, output_start]()->string {
          double*       typed_input  = static_cast<double*>(input_start);
          fftw_complex* typed_output = static_cast<fftw_complex*>(output_start);
          fftw_execute_dft_r2c(plan, typed_input, typed_output);
        };
      break;
    case FftProcessorConfig::FftType::kComplexToReal:
      thread_work = [plan, input_start, output_start]()->string {
          fftw_complex* typed_input  = static_cast<fftw_complex*>(input_start);
          double*       typed_output = static_cast<double*>(output_start);
          fftw_execute_dft_c2r(plan, typed_input, typed_output);
        };
      break;
    case FftProcessorConfig::FftType::kRealToReal:
      thread_work = [plan, input_start, output_start]()->string {
          double* typed_input  = static_cast<double*>(input_start);
          double* typed_output = static_cast<double*>(output_start);
          fftw_execute_r2r(plan, typed_input, typed_output);
        };
      break;
    default:
      cerr << "Unsupported FFT type '" << cfg_.type() << "'\n";
      abort();
    }

    // Sanity check; the default case above should have guaranteed this condition
    assert(thread_work);
    SubmitThreadWork(thread_work);
  }

  return "";
}
