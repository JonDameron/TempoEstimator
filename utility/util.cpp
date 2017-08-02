#include <cmath>
#include <string>
#include <unistd.h>
#include <sys/syscall.h>
#include "utility/util.hpp"

using namespace std;

void StackError (std::stringstream* err_strm, const std::string& str)
{
  if (err_strm->tellp() > 0) {
    err_strm->str(str + ": " + err_strm->str());
  } else {
    (*err_strm) << str;
  }
}

void AppendError (std::stringstream* err_strm, const std::string& str)
{
  if (err_strm->tellp() > 0) {
    (*err_strm) << "; " << str;
  } else {
    (*err_strm) << str;
  }
}

void AppendErrorWithErrno (std::stringstream* err_strm, const std::string& str)
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

void MakeHanningWindow (int n, std::vector<double>* out)
{
  if (out->size() < n) {
    out->resize(n);
  }

  // Save a little time by pre-calculating the constant scale factor in the
  // cos() argument
  const double cos_scale = 2.0 * M_PI / (n-1);

  for (int i = 0; i < n; ++i)
  {
    (*out)[i] = 0.5 * (1 - cos(i * cos_scale));
  }
}

pid_t GetCurrentThreadId ()
{
  return (pid_t) syscall(SYS_gettid);
}
