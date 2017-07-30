#pragma once

#include <memory>
#include <string>

#include "pipeline/abstract_pipeline_node.hpp"

namespace pipeline {

/** An implicit cast from FromType to ToType must be possible.
 * Since all class methods are templated, they're defined inline here instead
 * of in a separate cpp file to ensure that all needed implementations are
 * generated by the compiler.
 */
template <class FromType, class ToType>
class TypeConverter : public AbstractPipelineNode
{
public:

  TypeConverter ()
  {}

  ~TypeConverter ()
  {}

private:

  /** 1,000,000 type conversions per thread work unit is a miniscule amount
   * (a single thread on a modern CPU should be able to complete that many in a
   * small fraction of a second), but this is fine for demonstration purposes.
   */
  static const int kDefaultNWorkUnits = 1000*1000;

  std::string Process (std::shared_ptr<const ProcData> input,
                       std::shared_ptr<ProcData>* output) override
  {
    const size_t n_vals_total = input->size_bytes() / sizeof(FromType);
    auto in_buf = static_cast<const FromType*>(input->untyped_data());

    *output = ProcData::New(n_vals_total * sizeof(ToType));
    auto out_buf = static_cast<ToType*>((*output)->untyped_data());

    // [multithreading] The number of work units to be submitted for
    // parallel processing
    const int n_work_units = n_vals_total / kDefaultNWorkUnits;

    size_t val_offset = 0;

    for (int unit_i = 0; unit_i < n_work_units; ++unit_i)
    {
      // If n_vals_total doesn't evenly divide n_work_units, we'll have to
      // account for the remainder by distributing it across that many units.
      const int n_vals_this_unit = n_vals_total / n_work_units
          + (unit_i < n_vals_total % n_work_units ? 1 : 0);

      const FromType* in_buf_this_unit = in_buf + val_offset;

      ToType* out_buf_this_unit = out_buf + val_offset;

      auto thread_work_func =
          [n_vals_this_unit, in_buf_this_unit, out_buf_this_unit]()->std::string
          {
            for (size_t i = 0; i < n_vals_this_unit; ++i) {
              out_buf_this_unit[i] = static_cast<ToType>(in_buf_this_unit[i]);
            }
            return "";
          };

      SubmitThreadWork(thread_work_func);

      val_offset += n_vals_this_unit;
    }

    return "";
  }
};

}