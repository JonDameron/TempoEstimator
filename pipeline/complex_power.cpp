#include <complex>
#include <iostream>
#include <memory>
#include <string>
#include "pipeline/complex_power.hpp"
#include "utility/util.hpp"

using namespace std;
using namespace pipeline;

ComplexPower :: ComplexPower ()
{
}

ComplexPower :: ~ComplexPower ()
{
}

string ComplexPower :: Process (shared_ptr<const ProcData> input,
                                shared_ptr<ProcData>* output)
{
  stringstream res;

  const size_t n_vals_total = input->n_complex_doubles();
  const complex<double>* in_buf = input->complex_double_data();

  *output = ProcData::New(n_vals_total * sizeof(double));
  double* out_buf = (*output)->double_data();

  // The number of work units to be submitted for multithreaded processing
  const int n_work_units = kDefaultNWorkUnits;

  size_t val_offset = 0;

  for (int unit_i = 0; unit_i < n_work_units; ++unit_i)
  {
    // If n_vals_total doesn't evenly divide n_work_units, we'll have to
    // account for the remainder by distributing it across that many work units
    const int n_vals_this_unit = n_vals_total / n_work_units
        + (unit_i < n_vals_total % n_work_units ? 1 : 0);

    const complex<double>* in_buf_this_unit = in_buf + val_offset;

    double* out_buf_this_unit = out_buf + val_offset;

    auto thread_work_func =
        [n_vals_this_unit, in_buf_this_unit, out_buf_this_unit]
                                              (int thread_index)->std::string {
          (void)thread_index; // Avoid 'unused param' compiler warning
          for (size_t i = 0; i < n_vals_this_unit; ++i) {
            out_buf_this_unit[i] =
                Sqr(in_buf_this_unit[i].real()) + Sqr(in_buf_this_unit[i].imag());
          }
          return "";
        };

    SubmitThreadWork(thread_work_func);

    val_offset += n_vals_this_unit;
  }

  return "";
}
