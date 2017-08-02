#include "pipeline/abstract_pipeline_head.hpp"
#include "pipeline/audio_file_reader.hpp"

using namespace std;
using namespace pipeline;

AudioFileReader :: AudioFileReader (
    shared_ptr<ThreadSquad> thread_squad,
    shared_ptr<PcmMonoAudioData> existing_audio_data)
: AbstractPipelineHead(thread_squad),
  audio_data_(existing_audio_data)
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

  // Avoid copy overhead (uncompressed PCM audio files can be large) by invoking
  // std::move to transfer the single existing instance of the underlying data.
  *data_out = ProcData::New(audio_data_->MoveDataOut());

  // We no longer need the PCM audio file data pointed to by audio_data_.
  audio_data_.reset();

  return "";
}
