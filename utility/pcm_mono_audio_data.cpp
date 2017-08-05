#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include "utility/pcm_mono_audio_data.hpp"
#include "utility/util.hpp"

using namespace std;

const uint8_t PackedPcmAudioHeader::kFmt4cc  [4] = { 'f', 'm', 't', ' ' };
const uint8_t PackedPcmAudioHeader::kList4cc [4] = { 'L', 'I', 'S', 'T' };
const uint8_t PackedPcmAudioHeader::kRiff4cc [4] = { 'R', 'I', 'F', 'F' };
const uint8_t PackedPcmAudioHeader::kWave4cc [4] = { 'W', 'A', 'V', 'E' };
const uint8_t PackedPcmAudioHeader::kData4cc [4] = { 'd', 'a', 't', 'a' };

bool PackedPcmAudioHeader :: HasOptionalList () const
{
  return ( 0 == memcmp(
      this->optional_special_chunk_position, PackedPcmAudioHeader::kList4cc, 4) );
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
{
  // Read entire .wav file, then fill in header & data fields.

  // TODO: Perhaps delay our read of the data region, which could potentially be
  //     massive, until we actually need the data. The caller might only be
  //     interested in the header info.
  // TODO: Accept sample formats other than int16.

  string result;

  // Initialize packed struct class members
  memset(&packed_header_, 0, sizeof(packed_header_));

  // Open file in binary mode and immediately seek to end so we can easily
  // determine the total file size.

  ifstream file (wav_file_path, ios::binary | ios::ate);
  if (file.fail()) {
    AppendErrorWithErrno(&init_error_, "Failed to open audio file for reading");
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

  // The header format can vary significantly; it's a bit tricky to interpret
  // it properly

  const size_t initial_header_read_size =
      packed_header_.optional_special_chunk_content_position
      - reinterpret_cast<uint8_t*>(&packed_header_);

  file.read(reinterpret_cast<char*>(&packed_header_),
            initial_header_read_size);
  if (file.fail()) {
    AppendErrorWithErrno(&init_error_, "Failed to read audio file");
    return;
  }

  // Validate the header, and ensure the data format is one that we support.

  result = packed_header_.Validate();
  if (!result.empty()) {
    init_error_ << "Invalid RIFF/WAVE header: " << result;
    return;
  }

  result = CheckAudioFormatCompatibility(packed_header_);
  if (!result.empty()) {
    init_error_ << "Incompatible audio format: " << result;
    return;
  }

  const int n_bytes_per_sample = packed_header_.n_bits_per_sample / 8;

  // Some .wav headers ignore the cb_size, n_valid_bits_per_sample,
  // channel_mask, and subformat_info fields, and instead define a LIST chunk
  // at this position containing extra information, for example, the song name.

  if (packed_header_.HasOptionalList())
  {
    // Read LIST chunk, ignore remaining standard header fields
    const uint32_t list_size =
        *(uint32_t*)(packed_header_.optional_special_chunk_position + 4);

    // Calculate how much more of the audio file must be read to reach the
    // end of the LIST chunk
    const size_t remaining_header_size = list_size;

    ext_header_.resize(remaining_header_size);
    file.read(ext_header_.data(), remaining_header_size);
    if (file.fail()) {
      AppendErrorWithErrno(&init_error_,
          "Failed to read extra LIST data in file header");
      return;
    }
  }

  if (0 == memcmp(packed_header_.optional_special_chunk_position,
                  PackedPcmAudioHeader::kData4cc, 4))
  {
    memcpy( &packed_data_chunk_cap_,
            packed_header_.optional_special_chunk_position,
            sizeof(packed_data_chunk_cap_) );
  }
  else
  {
    // Read 'data' 4CC and 32-bit size
    file.read(reinterpret_cast<char*>(&packed_data_chunk_cap_),
              sizeof(packed_data_chunk_cap_));
    if (file.fail()) {
      AppendErrorWithErrno(&init_error_, "Failed to read 'data' chunk ID fields");
      return;
    }
    if (0 != memcmp(packed_data_chunk_cap_.id,
                    PackedPcmAudioHeader::kData4cc, 4)) {
      AppendErrorWithErrno(&init_error_, "Invalid 'data' chunk ID");
      return;
    }
  }

  // Note that the 'data' 4CC and 32-bit size are considered to be 'header' info
  const ios::pos_type total_header_size = file.tellg();

  const ios::pos_type remaining_file_size = total_file_size - total_header_size;

  // Read whatever audio samples are present, even if the actual amount of
  // data remaining in the file is less than the 'data' chunk size.

  const size_t interleaved_sample_unit_size =
      packed_header_.n_channels * n_bytes_per_sample;

  size_t data_chunk_size_bytes =
      std::min<size_t>(remaining_file_size, packed_data_chunk_cap_.size);

  data_chunk_size_bytes -= data_chunk_size_bytes % interleaved_sample_unit_size;

  const size_t n_total_samples_all_chans =
      data_chunk_size_bytes / n_bytes_per_sample;

  const int n_samples_per_chan =
      n_total_samples_all_chans / packed_header_.n_channels;

  const size_t chan_data_size_bytes =
      data_chunk_size_bytes / packed_header_.n_channels;

  sample_data_.reset(new SimdAlignedBuffer(chan_data_size_bytes));

  // If the audio is already mono (single channel), save time by simply reading
  // the data directly into our samples buffer
  if (1 == packed_header_.n_channels)
  {
    file.read(reinterpret_cast<char*>(sample_data_->data()),
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
    file.read(reinterpret_cast<char*>(orig_samples.data()),
              data_chunk_size_bytes);
    if (file.fail()) {
      AppendErrorWithErrno(&init_error_,
          "Failed to read data section of audio file");
      return;
    }

    ConvertMultiChannelSamplesToMono(
        packed_header_.n_channels, n_samples_per_chan, orig_samples.data(),
        static_cast<int16_t*>(sample_data_->data()) );
  }

  // <file> will be automatically closed on destruction
}

PcmMonoAudioData :: PcmMonoAudioData (const PackedPcmAudioHeader& header)
{
   packed_header_ = header;
}

PcmMonoAudioData :: PcmMonoAudioData (const PcmMonoAudioData& copy_src)
: packed_header_(copy_src.packed_header_),
  ext_header_(copy_src.ext_header_),
  packed_data_chunk_cap_(copy_src.packed_data_chunk_cap_),
  init_error_(copy_src.init_error_.str()),
  sample_data_(copy_src.sample_data_ ? new SimdAlignedBuffer(*copy_src.sample_data_)
                                     : nullptr)
{
}

PcmMonoAudioData :: ~PcmMonoAudioData ()
{
}

std::string PcmMonoAudioData :: SetSamples (
    int sample_offset, const std::vector<int16_t>& new_samples)
{
  if ((sample_offset + new_samples.size()) * sizeof(new_samples[0])
         > sample_data_->size()) {
    return "Assignment of new sample data would overflow the existing buffer";
  }

  memcpy( static_cast<int16_t*>(sample_data_->data()) + sample_offset,
          new_samples.data(),
          new_samples.size() * sizeof(new_samples[0]) );

  return "";
}

std::string PcmMonoAudioData :: OverdubSamples (
    int sample_offset, const std::vector<int16_t>& new_samples,
    double orig_audio_scaling)
{
  if ((sample_offset + new_samples.size()) * sizeof(new_samples[0])
         > sample_data_->size()) {
    return "Assignment of new sample data would overflow the existing buffer";
  }

  int16_t* samples = static_cast<int16_t*>(sample_data_->data());
  for (int i = 0; i < new_samples.size(); ++i) {
    samples[sample_offset + i] =
        samples[sample_offset + i] * orig_audio_scaling + new_samples[i];
  }

  return "";
}

string PcmMonoAudioData :: WriteToFile (const string& wav_file_path)
{
  ofstream file (wav_file_path, ios::out | ios::binary);
  if (file.fail()) {
    return "Failed to open audio file for writing";
  }

  PackedPcmAudioHeader out_header = packed_header_;
  out_header.n_channels = 1;
  out_header.riff_chunk_size -=
      sample_data_->size() * (packed_header_.n_channels-1);
  out_header.n_avg_bytes_per_sec /= packed_header_.n_channels;
  out_header.n_block_align /= packed_header_.n_channels;

  const size_t initial_header_write_size =
      packed_header_.optional_special_chunk_position
      - reinterpret_cast<uint8_t*>(&packed_header_);

  file.write(reinterpret_cast<const char*>(&out_header),
             initial_header_write_size);

  if (packed_header_.HasOptionalList())
  {
    file.write(
        reinterpret_cast<char*>(out_header.optional_special_chunk_position), 8);
    file.write(ext_header_.data(), ext_header_.size());
  }
  else if (0 != memcmp(out_header.optional_special_chunk_position,
                       PackedPcmAudioHeader::kData4cc, 4))
  {
    out_header.channel_mask = 1;
    out_header.n_valid_bits_per_sample =
        packed_header_.n_valid_bits_per_sample / packed_header_.n_channels;
    file.write(
        reinterpret_cast<char*>(out_header.optional_special_chunk_position),
        sizeof(packed_header_) - initial_header_write_size );
  }

  PackedDataChunkCap out_data_chunk_cap = packed_data_chunk_cap_;
  out_data_chunk_cap.size = sample_data_->size();
  file.write(reinterpret_cast<const char*>(&out_data_chunk_cap),
             sizeof(out_data_chunk_cap));

  file.write(reinterpret_cast<const char*>(sample_data_->data()),
             sample_data_->size());

  if (file.fail()) {
    return "Failed to write PCM audio header and sample data to file";
  }

  // <file> will be automatically closed on destruction

  return "";
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
