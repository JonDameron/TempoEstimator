#pragma once

#include <memory>
#include <string>

#include "pipeline/abstract_pipeline_head.hpp"
#include "utility/pcm_mono_audio_data.hpp"

namespace pipeline {

class AudioFileReader : public AbstractPipelineHead
{
public:

  /** Construct a new instance from existing PCM audio data.
   */
  explicit AudioFileReader (
      std::shared_ptr<ThreadSquad> thread_squad,
      std::shared_ptr<PcmMonoAudioData> existing_audio_data );

  // TODO: Add a constructor that takes a path to a PCM audio file and
  // attempts to initialize audio_data_ from the content therein.

  ~AudioFileReader ();

private:

  std::string CollectHeadData (std::shared_ptr<ProcData>* data_out) override;

  std::shared_ptr<PcmMonoAudioData> audio_data_;
};

}
