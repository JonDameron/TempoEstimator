#include "utility/proc_data.hpp"

using namespace std;

ProcData :: ProcData (size_t size_bytes)
: buf_(new SimdAlignedBuffer(size_bytes))
{
}

ProcData :: ProcData (unique_ptr<SimdAlignedBuffer>&& move_src)
: buf_(move_src)
{
}

ProcData :: ~ProcData ()
{
  // data_ is a unique_ptr, so the allocation will be automatically freed
}

shared_ptr<ProcData> ProcData :: New ()
{
  return make_shared<ProcData>(0);
}

shared_ptr<ProcData> ProcData :: New (size_t size_bytes)
{
  return make_shared<ProcData>(size_bytes);
}

shared_ptr<ProcData> ProcData :: New (unique_ptr<SimdAlignedBuffer>&& move_src)
{
  return make_shared<ProcData>(move_src);
}
