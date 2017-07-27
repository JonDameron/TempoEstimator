#pragma once

#include <complex>
#include <cstdint>
#include <cstdlib>
#include <memory>

class ProcData
{
private:

  /* New ProcData instances are expected to be constructed relatively
   * frequently. To avoid memory leaks, each instance must be wrapped in a
   * shared_ptr via the static method ProcData::New().
   */

  /** Constructs a new ProcData with uninitialized data buffer of given size.
   */
  explicit ProcData (size_t size_bytes);

  /** Transfers data from an existing SimdAlignedBuffer without incurring
   * copy overhead.
   */
  explicit ProcData (std::unique_ptr<SimdAlignedBuffer>&& move_src);

public:

  static std::shared_ptr<ProcData> New ();

  static std::shared_ptr<ProcData> New (size_t size_bytes);

  static std::shared_ptr<ProcData> New (std::unique_ptr<SimdAlignedBuffer>&& move_src);

  template <class ElementType>
  static std::shared_ptr<ProcData> New (size_t n_elements) {
    return std::make_shared<ProcData>(n_elements * sizeof(ElementType));
  }

  ~ProcData ();

  /** Convenience getter for buffer size that assumes an underlying type
   * of 'double'
   */
  size_t n_doubles () const {
    return buf_->size() / sizeof(double);
  }

  size_t size_bytes () const {
    return buf_->size();
  }

  const void* untyped_data () const {
    return buf_->data();
  }

  void* untyped_data () {
    return buf_->data();
  }

  const uint8_t* byte_data () const {
    return static_cast<const uint8_t*>(buf_->data());
  }

  uint8_t* byte_data () {
    return static_cast<uint8_t*>(buf_->data());
  }

  const double* double_data () const {
    return static_cast<const double*>(buf_->data());
  }

  double* double_data () {
    return static_cast<double*>(buf_->data());
  }

  const std::complex<double>* complex_double_data () const {
    return static_cast<const std::complex<double>*>(buf_->data());
  }

  std::complex<double>* complex_double_data () {
    return static_cast<std::complex<double>*>(buf_->data());
  }

private:

  /** By using unique_ptr here, we can std::move large data sets directly
   * to/from <data_> without incurring copy overhead.
   */
  std::unique_ptr<SimdAlignedBuffer> buf_;
};
