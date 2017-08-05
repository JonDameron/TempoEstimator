#pragma once

#include <cstring>
#include <sstream>
#include <vector>
#include <fftw3.h>

/** The override method qualifier keyword was unsupported on GCC prior to
 * version 4.7, even with the -std=c++0x option
 */
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
# define OVERRIDE override
#else
# define OVERRIDE
#endif

enum class DataType {
  kUndefined,
  kDouble,
  kComplexDouble
};

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
  : size_(size),
    data_(fftw_malloc(size))
  {}

  SimdAlignedBuffer (const SimdAlignedBuffer& copy_src)
  : size_(copy_src.size_),
    data_(memcpy(fftw_malloc(copy_src.size_), copy_src.data_, copy_src.size_))
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

template <class ValueType>
inline ValueType Sqr (ValueType val) {
  return val * val;
}

/** For concatenating a sequence of messages relating to a single error;
 * specifically intended to be used when unwinding a call stack on error.
 * @param err_strm Stream to receive the concatenated error messages.
 * @param str Error message to be added.
 * @returns Reference to *existing_err for convenience, as typically the result
 *    of this function is returned by the calling function
 */
std::string& StackError (std::string* existing_err, const std::string& new_err);

/** For concatenating error messages that aren't necessarily related.
 * @param err_strm Stream to receive the concatenated error messages.
 * @param str Error message to be added.
 */
void AppendError (std::stringstream* err_strm, const std::string& str);

/** Version of AppendError that includes strerror(errno) if errno is set.
 */
void AppendErrorWithErrno (std::stringstream* err_strm, const std::string& str);

/** Returns a unique ID for the invoking thread. Useful, e.g., for mappings
 * that need one unique entry per thread.
 * From gettid manual page:
 * gettid() returns the caller's thread ID (TID).  In a single-threaded
 * process, the thread ID is equal to the process ID (PID, as returned
 * by getpid(2)).  In a multithreaded process, all threads have the same
 * PID, but each one has a unique TID.
 * Glibc does not provide a wrapper for this system call; call it using
 * syscall(2).
 */
pid_t GetCurrentThreadId ();

bool IsSimdAligned (const void* addr);

/** Convert beats-per-minute to samples-per-beat
 */
double BpmToSamplesPerBeat (double bpm, double samples_per_sec);

/** Construct a uniform, a.k.a. rectangular window.
 * @param n Size of window.
 * @param out Output vector to receive window data.
 */
void MakeUniformWindow (int n, std::vector<double>* out);

/** Construct a standard Hanning window.
 * @param n Size of window.
 * @param out Output vector to receive window data.
 */
void MakeHanningWindow (int n, std::vector<double>* out);

/** Perform discrete unnormalized cross-correlation of two given input arrays.
 * @param corr_out_size Size, in # values (not bytes), of corr_out. This value
 * is used only to verify that the output array is large enough to store all
 * output data.
 */
std::string Correlate (const double* corr_in_1,
                       int corr_in_1_step,
                       int corr_in_1_len,
                       const double* corr_in_2,
                       int corr_in_2_step,
                       int corr_in_2_len,
                       double* corr_out,
                       int corr_out_len);
