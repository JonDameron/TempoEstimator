#pragma once

#include <sstream>
#include <vector>
#include <fftw3.h>


/** A simple wrapper class for fftw_malloc(). Faster SIMD
 * (Single Instruction Multiple Data) CPU instructions can often be utilized
 * for data that is specially aligned in memory -- typically on 16-byte
 * boundaries -- and so, especially on newer hardware, SIMD-aligning data
 * buffers can significantly boost performance. FFTW in particular takes
 * advantage of SIMD technology when possible.
 * Please see also:
 * http://www.fftw.org/fftw3_doc/SIMD-alignment-and-fftw_005fmalloc.html
 * Intel optimization manual: "Data must be 16-byte aligned when loading to and
 * storing from the 128-bit XMM registers used by SSE/SSE2/SSE3/SSSE3.
 * This must be done to avoid severe performance penalties."
 */
class SimdAlignedBuffer
{
public:

  explicit SimdAlignedBuffer (size_t size)
  : size_(size), data_(fftw_malloc(size))
  {}

  ~SimdAlignedBuffer ()
  { fftw_free(data_); }

  const void* data () const {
    return data_;
  }

  void* data () {
    return data_;
  }

  size_t size () const {
    return size_;
  }

private:

  size_t size_;
  void* data_;
};

/** For concatenating a sequence of messages relating to a single error;
 * specifcally intended to be used when unwinding a call stack on error.
 * @param err_strm Stream to receive the concatenated error messages.
 * @param str Error message to be added.
 */
void StackError (std::stringstream* err_strm, const std::string& str);

/** For concatenating error messages that aren't necessarily related.
 * @param err_strm Stream to receive the concatenated error messages.
 * @param str Error message to be added.
 */
void AppendError (std::stringstream* err_strm, const std::string& str);

/** Version of AppendError that includes strerror(errno) if errno is set.
 */
void AppendErrorWithErrno (std::stringstream* err_strm, const std::string& str);

/** Construct a standard Hanning window.
 * @param n Size of window.
 * @param out Output vector to receive window data.
 */
void MakeHanningWindow (int n, std::vector<double>* out);
