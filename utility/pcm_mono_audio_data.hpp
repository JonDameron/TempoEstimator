#pragma once

#include <cstdint>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>
#include "utility/util.hpp"

struct PackedPcmAudioHeader
{
  /** 4CC code for binary "RIFF" string at top of header
   */
  static const uint8_t kFmt4cc  [4];
  static const uint8_t kList4cc [4];
  static const uint8_t kRiff4cc [4];
  static const uint8_t kWave4cc [4];
  static const uint8_t kData4cc [4];

  bool HasOptionalList () const;

  std::string Validate () const;

  /** Must be "RIFF"
   */
  uint8_t riff_chunk_id [4];

  /** Does not account for sizes of chunk_id and chunk_size itself; thus,
   * chunk_size == total file size - 8.
   */
  uint32_t riff_chunk_size;

  /** Must be "WAVE"
   */
  uint8_t wave_id [4];

  /** Must be "fmt "
   */
  uint8_t fmt_chunk_id [4];

  /** Must be 16, 18, or 40
   */
  uint32_t chunk_size;

  uint16_t format_code;

  uint16_t n_channels;

  uint32_t n_samples_per_sec;

  uint32_t n_avg_bytes_per_sec;

  uint16_t n_block_align;

  uint16_t n_bits_per_sample;

  /* Some .wav headers ignore the cb_size, n_valid_bits_per_sample,
     channel_mask, and subformat_info fields, and instead define a LIST chunk
     at this position containing extra information, for example, the song name.
     This size-zero member acts as a pointer to the position of the header that
     might contain the LIST data instead of the 'standard' fields.
   */
  union
  {
    struct {
      /** .wav specification: "Size of the extension"
       */
      uint16_t cb_size;
      uint16_t n_valid_bits_per_sample;
      uint32_t channel_mask;
    } __attribute__((packed)) std_data;

    struct {
      uint8_t chunk_id [4];
      uint32_t chunk_size;
    } __attribute__((packed)) list_chunk_cap;

    struct {
      uint8_t chunk_id [4];
      uint32_t chunk_size;
    } __attribute__((packed)) data_chunk_cap;

    uint8_t bytes [8];
  }
  optional_special_chunk;

  uint8_t optional_special_chunk_content_position [0];

  /** Globally unique identifier (GUID), including the data format code
   */
  uint8_t subformat_info_ [16];
}
__attribute__((packed));


class PcmMonoAudioData
{
public:

  explicit PcmMonoAudioData (const std::string& wav_file_path);

  explicit PcmMonoAudioData (const PackedPcmAudioHeader& header);

  /** Copy constructor
   */
  PcmMonoAudioData (const PcmMonoAudioData& copy_src);

  ~PcmMonoAudioData ();

  std::string init_error () const {
    return init_error_.str();
  }

  const PackedPcmAudioHeader& header () const {
    return packed_header_;
  }

  size_t n_samples () const {
    return sample_data_->size() / (packed_header_.n_bits_per_sample / 8);
  }

  const int16_t* samples () const {
    return static_cast<const int16_t*>(sample_data_->data());
  }

  std::unique_ptr<SimdAlignedBuffer>&& MoveDataOut () {
    return std::move(sample_data_);
  }

  std::string SetSamples (int sample_offset,
                          const std::vector<int16_t>& new_samples);

  std::string OverdubSamples (int sample_offset,
                              const std::vector<int16_t>& new_samples,
                              double orig_audio_scaling = 1.0);

  std::string WriteToFile (const std::string& wav_file_path);

private:

  static std::string CheckAudioFormatCompatibility (
      const PackedPcmAudioHeader& header);

  /** Averages audio data across all interleaved channels of <in> to produce
   * mono output, then stores it in <mono_out>.
   * @param mono_out Output buffer for mono data. Size MUST be at least
   * sizeof(in[0]) * n_samples_per_chan.
   */
  static void ConvertMultiChannelSamplesToMono (
      int n_channels, size_t n_samples_per_chan, const int16_t* in,
      int16_t* mono_out );

  PackedPcmAudioHeader packed_header_;

  std::vector<char> ext_header_;

  struct PackedDataChunkCap {
    uint8_t id [4];
    uint32_t size;
  } __attribute((packed))
  packed_data_chunk_cap_;

  std::stringstream init_error_;

  std::unique_ptr<SimdAlignedBuffer> sample_data_;
};
