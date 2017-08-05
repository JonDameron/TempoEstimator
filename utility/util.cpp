#include <cmath>
#include <cstring>
#include <string>
#include <unistd.h>
#include <sys/syscall.h>
#include "utility/util.hpp"

using namespace std;

string& StackError (string* existing_err, const string& new_err)
{
  if (existing_err->empty() > 0) {
    *existing_err = new_err;
  } else {
    *existing_err = new_err + ": " + *existing_err;
  }
  return *existing_err;
}

void AppendError (stringstream* err_strm, const string& str)
{
  if (err_strm->tellp() > 0) {
    (*err_strm) << "; " << str;
  } else {
    (*err_strm) << str;
  }
}

void AppendErrorWithErrno (stringstream* err_strm, const string& str)
{
  // Invoke strerror(errno) before calling any function that may alter errno
  string errno_str = (errno ? strerror(errno) : "");

  if (str.empty()) {
    AppendError(err_strm, errno_str);
  } else if (errno_str.empty()) {
    AppendError(err_strm, str);
  } else {
    AppendError(err_strm, str + ": " + errno_str);
  }
}

pid_t GetCurrentThreadId ()
{
  return (pid_t) syscall(SYS_gettid);
}

bool IsSimdAligned (const void* addr)
{
  // TODO: 16 is the correct alignment constant for most systems, but not all.
  //       The literal should be replaced with a system-specific value.
  return  0 == reinterpret_cast<uint64_t>(addr) % 16;
}

double BpmToSamplesPerBeat (double bpm, double samples_per_sec)
{
  return samples_per_sec * 60.0 / bpm;
}

void MakeUniformWindow (int n, vector<double>* out)
{
  // Note that vector::clear() does NOT alter capacity, so we aren't risking
  // an unnecessary deallocation / reallocation
  out->clear();
  out->resize(n, 1.0);
}

void MakeHanningWindow (int n, vector<double>* out)
{
  if (out->size() < n) {
    out->resize(n);
  }

  // Save a little time by pre-calculating the constant scale factor in the
  // cos() argument
  const double cos_scale = 2.0 * M_PI / (n-1);

  for (int i = 0; i < n; ++i) {
    (*out)[i] = 0.5 * (1.0 - cos(i * cos_scale));
  }
}

string Correlate (const double* corr_in_1,
                  int corr_in_1_step,
                  int corr_in_1_len,
                  const double* corr_in_2,
                  int corr_in_2_step,
                  int corr_in_2_len,
                  double* corr_out,
                  int corr_out_len)
{
  stringstream strm;

  if (corr_in_1_len < corr_in_2_len) {
    strm << "Input-1 length '" << corr_in_1_len
         << "' must be >= input-2 length '" << corr_in_2_len << "'";
    return strm.str();
  }

  const int comp_len = corr_in_1_len - corr_in_2_len + 1;

  if (corr_out_len < comp_len) {
    strm << "The length '" << corr_out_len
         << "' of the output buffer is insufficient to accommodate the total"
         << " output length '" << comp_len << "' of the correlation";
    return strm.str();
  }

  for (int out_i = 0; out_i < comp_len; ++out_i)
  {
    double sum = 0;
    for (int in2_i = 0; in2_i < corr_in_2_len; ++in2_i)
    {
      // Note that out_i is also our current input-1 offset
      sum +=   corr_in_1[(out_i + in2_i) * corr_in_1_step]
             * corr_in_2[in2_i * corr_in_2_step];
    }
    corr_out[out_i] = sum;
  }

  return "";
}
