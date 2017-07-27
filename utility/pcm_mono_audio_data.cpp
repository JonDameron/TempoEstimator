#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include "utility/pcm_mono_audio_data.hpp"
#include "utility/util.hpp"

using namespace std;

bool PackedPcmAudioHeader :: HasOptionalList () const
{
  return ( 0 == memcmp(
      this->optional_list_chunk_position, PackedPcmAudioHeader::kList4cc, 4) );
}

string PackedPcmAudioHeader :: Validate () const
{
  stringstream err;

  if (0 != memcmp(&this->riff_chunk_id, kRiff4cc, 4)) {
    AppendError(&err, "Invalid RIFF chunk ID");
  }
  if (0 != memcmp(&this->wave_id, kWave4cc, 4)) {
    AppendError(&err, "Invalid WAVE ID");
  }
  if (0 != memcmp(&this->fmt_chunk_id, kFmt4cc, 4)) {
    AppendError(&err, "Invalid 'fmt ' chunk ID");
  }
  // TODO: Perform additional validation

  return err.str();
}


PcmMonoAudioData :: PcmMonoAudioData (const string& wav_file_path)
: PcmMonoAudioData(wav_file_path, 1)
{
}

PcmMonoAudioData :: PcmMonoAudioData (const string& wav_file_path,
                                      int n_samples_align)
{
  // Read entire .wav file, then fill in header & data fields.

  // TODO: Perhaps delay our read of the data region, which could potentially be
  //     massive, until we actually need the data. The caller might only be
  //     interested in the header info.
  // TODO: Accept sample formats other than int16.

  string result;

  // Open file in binary mode and immediately seek to end so we can easily
  // determine the total file size.

  ifstream file (wav_file_path, ios::binary | ios::ate);
  if (file.fail()) {
    AppendErrorWithErrno(&init_error_, "Failed to open audio file");
    return;
  }

  const ios::pos_type total_file_size = file.tellg();

  if (total_file_size < sizeof(PackedPcmAudioHeader)) {
    init_error_ << "Audio file is too small to contain a RIFF/WAVE header (size = "
        << total_file_size << ")";
    return;
  }

  // Seek to beginning of file, then read the whole thing.

  file.seekg(0, ios::beg);

  file.read(reinterpret_cast<ios::char_type*>(&header_), sizeof(header_));
  if (file.fail()) {
    AppendErrorWithErrno(&init_error_, "Failed to read audio file");
    return;
  }

  // Validate the header, and ensure the data format is one that we support.

  result = header_.Validate();
  if (!result.empty()) {
    init_error_ << "Invalid RIFF/WAVE header: " << result;
    return;
  }

  result = CheckAudioFormatCompatibility(header_);
  if (!result.empty()) {
    init_error_ << "Incompatible audio format: " << result;
    return;
  }

  const int n_bytes_per_sample = header_.n_bits_per_sample / 8;

  // Some .wav headers ignore the cb_size, n_valid_bits_per_sample,
  // channel_mask, and subformat_info fields, and instead define a LIST chunk
  // at this position containing extra information, for example, the song name.

  if (header_.HasOptionalList())
  {
    // Read LIST chunk, ignore remaining standard header fields
    const uint32_t list_size =
        *(uint32_t*)(header_.optional_list_chunk_position + 4);

    // Subtract an additional 8 bytes to account for size of LIST 4CC and
    // 32-bit LIST size that immediately follows
    const size_t list_data_already_read_size = sizeof(PackedPcmAudioHeader)
        - (header_.optional_list_chunk_position - &header_) - 8;

    // Calculate how much more of the audio file must be read to reach the
    // end of the LIST chunk
    const size_t remaining_header_size =
        list_size - list_data_already_read_size;

    ext_header_.resize(remaining_header_size);
    file.read(ext_header_.data(), remaining_header_size);
    if (file.fail()) {
      AppendErrorWithErrno(&init_error_,
          "Failed to read extra LIST data in file header");
      return;
    }
  }

  // Read 'data' 4CC and immediately following 32-bit size
  file.read(reinterpret_cast<ios::char_type*>(&data_chunk_cap_),
            sizeof(data_chunk_cap_));
  if (file.fail()) {
    AppendErrorWithErrno(&init_error_,
        "Failed to read 'data' chunk ID fields");
    return;
  }

  // Note that the 'data' 4CC and 32-bit size are considered to be 'header' info
  const ios::pos_type total_header_size = file.tellg();

  const ios::pos_type remaining_file_size = total_file_size - total_header_size;

  // Read whatever audio samples are present, even if the actual amount of
  // data remaining in the file is less than the 'data' chunk size.

  const size_t data_chunk_size_bytes =
      std::min<size_t>(remaining_file_size, data_chunk_cap_.size);

  const size_t interleaved_sample_unit_size =
      header_.n_channels * n_bytes_per_sample;

  data_chunk_size_bytes -= data_chunk_size_bytes % interleaved_sample_unit_size;

  const size_t n_total_samples_all_chans =
      data_chunk_size_bytes / n_bytes_per_sample;

  const int n_samples_per_chan = n_total_samples_all_chans / header_.n_channels;

  const size_t chan_data_size_bytes = data_chunk_size_bytes / header_.n_channels;

  n_samples_ = n_samples_per_chan;

  sample_data_.reset(new SimdAlignedBuffer(chan_data_size_bytes));

  // If the audio is already mono (single channel), save time by simply reading
  // the data directly into our samples buffer
  if (1 == header_.n_channels)
  {
    file.read(reinterpret_cast<ios::char_type*>(sample_data_.get()),
              data_chunk_size_bytes);
    if (file.fail()) {
      AppendErrorWithErrno(&init_error_,
          "Failed to read data section of mono audio file");
      return;
    }
  }
  else
  {
    vector<int16_t> orig_samples (n_total_samples_all_chans);
    file.read(reinterpret_cast<ios::char_type*>(orig_samples.data()),
              data_chunk_size_bytes);
    if (file.fail()) {
      AppendErrorWithErrno(&init_error_,
          "Failed to read data section of audio file");
      return;
    }

    ConvertMultiChannelSamplesToMono(
        header_.n_channels, n_samples_per_chan, orig_samples.data(),
        static_cast<int16_t*>(sample_data_.get()) );
  }
}

PcmMonoAudioData :: ~PcmMonoAudioData ()
{
}

string PcmMonoAudioData :: CheckAudioFormatCompatibility (
    const PackedPcmAudioHeader& header)
{
  if (16 != header.n_bits_per_sample) {
    return "n_bits_per_sample is '" + to_string(header.n_bits_per_sample)
        + "', must be 16";
  }
  if (header.n_channels > 2) {
    return "n_channels is '" + to_string(header.n_channels)
        + "', must be <= 2";
  }

  return "";
}

void PcmMonoAudioData :: ConvertMultiChannelSamplesToMono (
    int n_channels, size_t n_samples_per_chan, const int16_t* in,
    int16_t* mono_out)
{
  if (1 == n_channels) {
    // Already mono; caller should probably handle this case with a std::move
    // if possible, especially if <n_samples_per_chan> is very large
    memcpy(mono_out, in, n_samples_per_chan * sizeof(in[0]));
    return;
  }

  // TODO: Possibly support n_channels > 2 someday?
  assert(2 == n_channels);

  for (int mono_i = 0; mono_i < n_samples_per_chan; ++mono_i)
  {
    const int multi_i = mono_i * n_channels;
    mono_out[mono_i] = (in[multi_i] + in[1+multi_i]) / 2;
  }
}
