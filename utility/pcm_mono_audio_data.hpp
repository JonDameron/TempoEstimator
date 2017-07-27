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
  static const uint8_t kFmt4cc  [4] = { 'f', 'm', 't', ' ' };
  static const uint8_t kList4cc [4] = { 'L', 'I', 'S', 'T' };
  static const uint8_t kRiff4cc [4] = { 'R', 'I', 'F', 'F' };
  static const uint8_t kWave4cc [4] = { 'W', 'A', 'V', 'E' };

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
     Since a LIST at this position can be any size (specified after the 'LIST'
     4CC), we can't make this cleaner via a union.
   */
  uint8_t optional_list_chunk_position [0];

  /** .wav specification: "Size of the extension"
   */
  uint16_t cb_size;

  uint16_t n_valid_bits_per_sample;

  uint32_t channel_mask;

  /** Globally unique identifier (GUID), including the data format code
   */
  uint8_t subformat_info_ [16];
}
__attribute__((packed));


class PcmMonoAudioData
{
public:

  explicit PcmMonoAudioData (const std::string& wav_file_path);

  PcmMonoAudioData (const std::string& wav_file_path, int n_samples_align);

  ~PcmMonoAudioData ();

  const std::string& init_error () const {
    return init_error_.str();
  }

  const PackedPcmAudioHeader& header () const {
    return header_;
  }

  const size_t n_samples () const {
    return n_samples_;
  }

  const int16_t* samples () const {
    return (const int16_t*)sample_data_.get();
  }

  std::unique_ptr<SimdAlignedBuffer>&& MoveDataOut () {
    return std::move(sample_data_);
  }

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

  std::stringstream init_error_;

  PackedPcmAudioHeader header_;

  std::vector<std::ios::char_type> ext_header_;

  struct { uint8_t id [4]; uint32_t size; } __attribute__((packed)) data_chunk_cap_;

  size_t n_samples_;

  std::unique_ptr<SimdAlignedBuffer> sample_data_;
};
