#include <fftw3.h>
#include <string>
#include "utility/util.hpp"
#include "pipeline/fft_processor_config.hpp"

using namespace std;
using namespace pipeline;

string pipeline::FftTypeToString (FftType type)
{
  // Using simple 'c2r'-style to follow convention established by FFTW API names
  switch (type) {
  case FftType::kComplexToComplex: return "c2c";
  case FftType::kComplexToReal:    return "c2r";
  case FftType::kRealToComplex:    return "r2c";
  case FftType::kRealToReal:       return "r2r";
  default: break;
  }
  return "unknown";
}


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
    AppendError(&err, "type '" + to_string(static_cast<int>(type_))
        + "' is invalid");
  }
  if (fft_len_ <= 0) {
    AppendError(&err, "fft_len '" + to_string(fft_len_) + "' must be > 0");
  }
  if (n_ffts_ <= 0) {
    AppendError(&err, "n_ffts '" + to_string(n_ffts_) + "' must be > 0");
  }
  if (in_stride_ <= 0) {
    AppendError(&err, "in_stride '" + to_string(in_stride_) + "' must be > 0");
  }
  if (in_dist_ <= 0) {
    AppendError(&err, "in_dist '" + to_string(in_dist_) + "' must be > 0");
  }
  if (out_stride_ <= 0) {
    AppendError(&err, "out_stride '" + to_string(out_stride_) + "' must be > 0");
  }
  if (out_dist_ <= 0) {
    AppendError(&err, "out_dist '" + to_string(out_dist_) + "' must be > 0");
  }
  if (FFTW_FORWARD != sign_ && FFTW_BACKWARD != sign_) {
    AppendError(&err, "sign '" + to_string(sign_)
        + "' must be FFTW_FORWARD or FFTW_BACKWARD");
  }
  if (fftw_planner_time_limit_sec_ < 0
        && fftw_planner_time_limit_sec_ != FFTW_NO_TIMELIMIT) {
    AppendError(&err, "fftw_planner_time_limit_sec '"
        + to_string(fftw_planner_time_limit_sec_)
        + "' must be >= 0 or FFTW_NO_TIMELIMIT");
  }

  return err.str();
}

string FftProcessorConfig :: ToString () const
{
  stringstream strm;

  strm << "type = " << FftTypeToString(type_) << "\n";
  strm << "fft_len = " << fft_len_ << "\n";
  strm << "n_ffts = " << n_ffts_ << "\n";
  strm << "in_stride = " << in_stride_ << "\n";
  strm << "in_dist = " << in_dist_ << "\n";
  strm << "out_stride = " << out_stride_ << "\n";
  strm << "out_dist = " << out_dist_ << "\n";
  strm << "sign = " << sign_ << "\n";
  strm << "fftw_planner_flags = " << fftw_planner_flags_ << "\n";
  //strm << "fftw_planner_time_limit_sec = " << fftw_planner_time_limit_sec_ << "\n";
  strm << "real_to_real_kind = " << real_to_real_kind_ << "\n";
  strm << "accept_odd_fft_len_for_real_transform = "
       << (accept_odd_fft_len_for_real_transform_ ? "true" : "false");

  return strm.str();
}
