#include "pipeline/abstract_pipeline_head.hpp"
#include "pipeline/audio_file_reader.hpp"

using namespace std;
using namespace pipeline;

AudioFileReader :: AudioFileReader (
    shared_ptr<ThreadSquad> thread_squad,
    unique_ptr<PcmMonoAudioData>&& existing_audio_data)
: AbstractPipelineHead(thread_squad),
  audio_data_(std::move(existing_audio_data))
{
}

AudioFileReader :: ~AudioFileReader ()
{
}

string AudioFileReader :: CollectHeadData (shared_ptr<ProcData>* data_out)
{
  if (!audio_data_) {
    // This method has already been invoked, and this particular
    // AbstractPipelineHead implementation only returns non-empty data once
    return "";
  }
  if (!audio_data_->init_error().empty()) {
    return "Cannot collect audio data, there are initialization error(s): "
        + audio_data_->init_error();
  }

  //*data_out = ProcData::New(audio_data_->MoveDataOut());
  const size_t audio_data_size =
      sizeof(audio_data_->samples()[0]) * audio_data_->n_samples();
  *data_out = ProcData::New(audio_data_size);
  memcpy((*data_out)->untyped_data(), audio_data_->samples(), audio_data_size);

  return "";
}
