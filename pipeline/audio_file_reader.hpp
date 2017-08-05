#pragma once

#include <memory>
#include <string>

#include "pipeline/abstract_pipeline_head.hpp"
#include "utility/pcm_mono_audio_data.hpp"
#include "utility/util.hpp"

namespace pipeline {

class AudioFileReader : public AbstractPipelineHead
{
public:

  /** Construct a new instance from existing PCM audio data.
   * @param thread_squad Existing ThreadSquad instance to receive thread work
   * units associated with pipeline data processing.
   * @existing_audio_data Existing PCM mono audio data. Note that this is
   * required to be wrapped in a unique_ptr, the managed object of which will be
   * moved into the AudioFileReader. This is to avoid unnecessarily copying
   * large amounts of raw audio data.
   * After the constructor returns, the PCM audio header info is still
   * accessible via the audio_header() getter.
   */
  AudioFileReader (
      std::shared_ptr<ThreadSquad> thread_squad,
      std::unique_ptr<PcmMonoAudioData>&& existing_audio_data );

  // TODO: Add a constructor that takes a path to a PCM audio file and
  // attempts to initialize audio_data_ from the content therein.

  ~AudioFileReader ();

  const PackedPcmAudioHeader& audio_header () const {
    return audio_data_->header();
  }

  std::shared_ptr<PcmMonoAudioData> ClonePcmAudioData () {
    return std::make_shared<PcmMonoAudioData>(*audio_data_);
  }

private:

  std::string CollectHeadData (std::shared_ptr<ProcData>* data_out) OVERRIDE;

  std::unique_ptr<PcmMonoAudioData> audio_data_;
};

}
