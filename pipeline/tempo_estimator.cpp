#include <string>
#include <memory>
#include "pipeline/tempo_estimator.hpp"

using namespace std;
using namespace pipeline;

TempoEstimator :: TempoEstimator ()
{
}

TempoEstimator :: ~TempoEstimator ()
{
}

string TempoEstimator :: Process (shared_ptr<const ProcData> input,
                                  shared_ptr<ProcData>* output)
{
  // Instead of setting *output, TempoEstimator sets internal variables.
  // The final calculations are accessible via public getters.
  // The caller will interpret *output being unset as an indication that this
  // derived pipeline node is an endpoint.

  // <input> must be the real-valued output of the second-level FFT.



  return "";
}
