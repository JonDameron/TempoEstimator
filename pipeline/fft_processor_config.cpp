#include <fftw3.h>
#include "utility/util.hpp"
#include "pipeline/fft_processor_config.hpp"

using namespace std;
using namespace pipeline;

FftProcessorConfig :: FftProcessorConfig ()
: type_(FftType::kInvalid),
  fft_len_(0),
  n_ffts_(0),
  in_stride_(0),
  in_dist_(0),
  out_stride_(0),
  out_dist_(0),
  sign_(FFTW_FORWARD),
  fftw_planner_flags_(FFTW_MEASURE),
  fftw_planner_time_limit_sec_(kInvalidFftwPlannerTimeLimit),
  real_to_real_kind_(FFTW_DHT), // discrete Hartley r2r transform by default
  accept_odd_fft_len_for_real_transform_(false)
{
}

FftProcessorConfig :: ~FftProcessorConfig ()
{
}

string FftProcessorConfig :: Validate () const
{
  stringstream err;

  if (type_ < FftType::kMin || type_ > FftType::kMax) {
    AppendErr(&err, "type '" + to_string(static_cast<int>(type_))
        + "' is invalid");
  }
  if (fft_len_ <= 0) {
    AppendErr(&err, "fft_len '" + to_string(fft_len_) + "' must be > 0");
  }
  if (n_ffts_ <= 0) {
    AppendErr(&err, "n_ffts '" + to_string(n_ffts_) + "' must be > 0");
  }
  if (in_stride_ <= 0) {
    AppendErr(&err, "in_stride '" + to_string(in_stride_) + "' must be > 0");
  }
  if (in_dist_ <= 0) {
    AppendErr(&err, "in_dist '" + to_string(in_dist_) + "' must be > 0");
  }
  if (out_stride_ <= 0) {
    AppendErr(&err, "out_stride '" + to_string(out_stride_) + "' must be > 0");
  }
  if (out_dist_ <= 0) {
    AppendErr(&err, "out_dist '" + to_string(out_dist_) + "' must be > 0");
  }
  if (FFTW_FORWARD != sign_ && FFTW_BACKWARD != sign_) {
    AppendErr(&err, "sign '" + to_string(sign)
        + "' must be FFTW_FORWARD or FFTW_BACKWARD");
  }
  if (fftw_planner_time_limit_sec_ < 0
        && fftw_planner_time_limit_sec_ != FFTW_NO_TIMELIMIT) {
    AppendErr(&err, "fftw_planner_time_limit_sec '"
        + to_string(fftw_planner_time_limit_sec_)
        + "' must be >= 0 or FFTW_NO_TIMELIMIT");
  }

  return err.str();
}
