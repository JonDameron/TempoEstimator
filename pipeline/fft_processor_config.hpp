#pragma once

#include <vector>

namespace pipeline {

class FftProcessorConfig
{
public:

  enum class FftType {
    kInvalid = 0,
    kMin,
    kComplexToComplex = kMin,
    kComplexToReal,
    kRealToComplex,
    kRealToReal,
    kMax = kRealToReal
  };

  /** Special time limit value that prevents fftw_set_timelimit() from being
   * invoked at all. It's important to support a means of disabling this
   * functionality since FFTW uses a global setting for the planner time limit,
   * so multiple instances of FftProcessor using different time limit settings
   * may override one another. The application must be careful to set this
   * for just one instance of FftProcessor.
   * The FFTW manual states that special valid time limit FFTW_NO_TIMELIMIT
   * is negative, so subtracting 1 to obtain another special value not equal
   * to FFTW_NO_TIMELIMIT.
   */
  static constexpr double kInvalidFftwPlannerTimeLimit = FFTW_NO_TIMELIMIT - 1.0;

  FftProcessorConfig ();

  ~FftProcessorConfig ();

  /** Type of DFT being performed
   */
  FftType type () const {
    return type_;
  }
  void set_type (FftType val) {
    type_ = val;
  }

  /** Length of each FFT
   */
  int fft_len () const {
    return fft_len_;
  }
  void set_fft_len (int val) {
    fft_len_ = val;
  }

  /** Number of FFTs to execute
   */
  int n_ffts () const {
    return n_ffts_;
  }
  void set_n_ffts (int val) {
    n_ffts_ = val;
  }

  /** Distance from one time-domain input sample to the next
   */
  int in_stride () const {
    return in_stride_;
  }
  void set_in_stride (int val) {
    in_stride_ = val;
  }

  /** Distance between consecutive FFT input sets; that is, from time-domain
   * input sample[N] of FFT[M] to sample[N] of FFT[M+1].
   */
  int in_dist () const {
    return in_dist_;
  }
  void set_in_dist (int val) {
    in_dist_ = val;
  }

  /** Distance from one frequency-domain output bin to the next
   */
  int out_stride () const {
    return out_stride_;
  }
  void set_out_stride (int val) {
    out_stride_ = val;
  }

  /** Distance between consecutive FFT output sets; that is, from
   * frequency-domain output bin[N] of FFT[M] to bin[N] of FFT[M+1]
   */
  int out_dist () const {
    return out_dist_;
  }
  void set_out_dist (int val) {
    out_dist_ = val;
  }

  int fftw_planner_flags () const {
    return fftw_planner_flags_;
  }
  void set_fftw_planner_flags (int val) {
    fftw_planner_flags_ = val;
  }

  double fftw_planner_time_limit_sec () const {
    return fftw_planner_time_limit_sec_;
  }

  /** FFTW uses a global setting for the planner time limit, so multiple
   * instances of FftProcessor using different time limit settings may override
   * one another. The application must be careful to invoke this function with
   * a value != kInvalidFftwPlannerTimeLimit for only a single FftProcessor.
   */
  void set_fftw_planner_time_limit_sec (double val) {
    fftw_planner_time_limit_sec_ = val;
  }

  /** Taper weights to apply to input data prior to FFT evaluation
   */
  const std::vector<double>& window () const {
    return window_;
  }
  void set_window (const std::vector<double>& val) {
    window_ = val;
  }

  int sign () const {
    return sign_;
  }
  void set_sign (int val) {
    sign_ = val;
  }

  int total_input_span_n_values () const {
    return in_dist_ * (n_ffts_-1) + in_stride_ * (fft_len_-1) + 1;
  }

  int total_output_span_n_values () const {
    return out_dist_ * (n_ffts_-1) + out_stride_ * (fft_len_-1) + 1;
  }

  size_t in_value_size () const {
    return type_ == FftType::kComplexToComplex || type_ == FftType::kComplexToReal
           ? sizeof(fftw_complex)
           : sizeof(double);
  }

  size_t out_value_size () const {
    return type_ == FftType::kComplexToComplex || type_ == FftType::kRealToComplex
           ? sizeof(fftw_complex)
           : sizeof(double);
  }

  /** If false (the default), the FftProcessor will reject a configuration that
   * would transform to or from real-valued input AND has an odd fft_len
   * AND would span more than one thread work unit.  This is because we want
   * to maximize efficiency by taking advantage of SIMD-aligned data, which
   * requires (on most machines) the input and output buffers to be allocated
   * on a 16-byte boundary.  To split a single contiguous in/out buffer up for
   * processing by multiple work units without compromising the 16-byte
   * boundary condition,  the underlying type must be fftw_complex (2 doubles)
   * OR, for real-valued data (1 double per datum), fft_len MUST be even.
   */
  void set_accept_odd_fft_len_for_real_transform (bool val) {
    accept_odd_fft_len_for_real_transform_ = val;
  }

  std::string Validate () const;

private:

  // See documentation for getters above

  FftType type_;
  int fft_len_;
  int n_ffts_;
  int in_stride_;
  int in_dist_;
  int out_stride_;
  int out_dist_;
  int sign_;
  int fftw_planner_flags_;
  double fftw_planner_time_limit_sec_;
  bool accept_odd_fft_len_for_real_transform_;

  std::vector<double> window_;
};

}
