#include <cstring>
#include "utility/proc_data.hpp"

using namespace std;

ProcData :: ProcData (size_t size_bytes)
: buf_(new SimdAlignedBuffer(size_bytes))
{
}

ProcData :: ProcData (const ProcData& copy_src)
: buf_(new SimdAlignedBuffer(copy_src.size_bytes()))
//: ProcData(copy_src.size_bytes()) older GCC doesn't support ctor delegation
{
  memcpy(buf_->data(), copy_src.untyped_data(), copy_src.size_bytes());
}

ProcData :: ProcData (unique_ptr<SimdAlignedBuffer>&& move_src)
: buf_(std::move(move_src))
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

shared_ptr<ProcData> ProcData :: New (const ProcData& copy_src)
{
  return make_shared<ProcData>(copy_src);
}

shared_ptr<ProcData> ProcData :: New (unique_ptr<SimdAlignedBuffer>&& move_src)
{
  return make_shared<ProcData>(std::move(move_src));
}
