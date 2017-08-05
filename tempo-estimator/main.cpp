/* Musical Tempo Estimation Utility
 * (c) 2017 Jonathan G. Dameron | jon.g.dameron@gmail.com
 */

#include <climits>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "utility/pcm_mono_audio_data.hpp"
#include "utility/util.hpp"
#include "pipeline/audio_file_reader.hpp"
#include "pipeline/complex_power.hpp"
#include "pipeline/fft_processor.hpp"
#include "pipeline/fft_processor_config.hpp"
#include "pipeline/tempo_estimator.hpp"
#include "pipeline/type_converter.hpp"

using namespace std;
using pipeline::FftType;

/*
 * TODO: Explore other tempo estimation approaches.
 * For example, another way to accomplish this is via autocorrelation of the
 * input audio to identify repetitions of a musical pattern. There are several
 * ways that computational requirements could be reduced, among them sample
 * decimation (process only every Nth audio sample) and autocorrelation of only
 * a small portion of the input.
 *
 * A major advantage of the current method, which performs a series of
 * time-domain FFTs on the input followed by a series of corner-turn FFTs on
 * the output of the first FFT, is that it primarily relies on the FFTW
 * library, which is heavily optimized to perform exceptionally well on a wide
 * variety of architectures.
 */

/** Default number of FFTs to execute per second while processing the raw
 * audio data.  Setting this too low will reduce the accuracy of the tempo
 * calculation; 200 < value < (sampleRate / 50) will probably generate a
 * reasonably accurate estimate.
 */
static const double kDefaultLevelOneFftsPerSecond = 200;

/** Since the level-2 FFT is a cornerturn (consecutive input data aren't
 * contiguous in memory), it isn't essential that kDefaultLevelTwoFftLen be
 * an even number to guarantee future in/out buffer SIMD alignment, unlike the
 * level-1 FFT length and kDefaultLevelTwoFftSpacing.
 * This duration should be long enough to include most any off-tempo song intro
 * and a substantial amount of regular-tempo music.
 */
static const double kDefaultLevelTwoFftSpanSeconds = 30;

/** MUST be an even number to guarantee future SIMD alignment of in/out buffers.
 */
static const int kDefaultLevelTwoFftSpacing = 2; //16

static const int kDefaultNPipelineProcThreads = 4;

static const double kDefaultFftwPlannerTimeLimitSec = 1.0;

static const double kMinAcceptableTempoBpm = 61.0;
static const double kMaxAcceptableTempoBpm = 181.0;

static const double kMetronomeTickAmplitude   = 0.30;
static const double kMetronomeTickDurationSec = 0.12;
static const double kMetronomeToneFreq        = 220; // 220 => A lovely 'A3'

/** Sadly, because there is no simple way that I am aware of to pass arbitrary
 * user data to a signal handler (not even via sigaction), we must use a global
 * as a workaround.
 */
shared_ptr<pipeline::AbstractPipelineHead> sighandler_data_pipeline_head;


static void Usage (const string& progname)
{
  cout << "Music file tempo estimation utility\n"
       << "(c) 2017 Jonathan G. Dameron\n"
       << "USAGE: " << progname << "\n"
       << "  -f  Set frequency profile granularity (N per second) [default = "
           << kDefaultLevelOneFftsPerSecond << "]\n"
       << "  -i  Set raw PCM input file (.wav)\n"
       << "  -h  Show this usage message\n"
       ;
}

static void HandleSignal (int signo)
{
  static int termsig_count = 0;

  if (SIGTERM == signo || SIGINT == signo)
  {
    ++termsig_count;
    if (1 == termsig_count)
    {
      cerr << "Caught signal '" << strsignal(signo)
           << "', cleaning up and quitting\n";
      // Attempt to acquire a local reference to the managed object before
      // performing any operations on it, since the global reference may be
      // nullified at any time in a multithreaded environment.
      shared_ptr<pipeline::AbstractPipelineHead> pipeline_head_local_ref =
          sighandler_data_pipeline_head;
      if (pipeline_head_local_ref) {
        pipeline_head_local_ref->Stop();
      }
    }
    else if (2 == termsig_count)
    {
      cerr << "Caught multiple termination signals, exiting\n";
      // Delay exit() for just a second to give any remaining threads a chance
      // to terminate cleanly; e.g., perhaps the user simply pressed Ctrl+C
      // a few times in rapid succession.
      sleep(1);
      exit(0);
    }
    // Ignore further termination signals (we've already handled the first two)
  }
  else
  {
    cerr << "Caught unhandled signal '" << strsignal(signo) << "'\n";
  }
}

void OverdubMetronome(
    double tempo_offset_seconds,
    double tempo_bpm,
    double tick_amplitude,
    double tone_frequency,
    double tick_duration_sec,
    PcmMonoAudioData* out_audio )
{
  // TODO: Eventually support sample formats other than signed 16-bit int
  assert(16 == out_audio->header().n_bits_per_sample);

  const int metronome_n_samples =
      (int)(tick_duration_sec * out_audio->header().n_samples_per_sec);

  vector<int16_t> metronome_tick (metronome_n_samples);
  for (int i = 0; i < metronome_n_samples; ++i) {
    metronome_tick[i] =
        tick_amplitude * SHRT_MAX
        * cos(2*M_PI*tone_frequency * (double)i / metronome_n_samples);
  }

  const double tempo_n_samples_per_beat =
      BpmToSamplesPerBeat(tempo_bpm, out_audio->header().n_samples_per_sec);

  // -1 to avoid slightly overrunning the output buffer
  const int n_ticks =
      (int)(out_audio->n_samples() / tempo_n_samples_per_beat) - 1;

  const int tempo_offset_n_samples =
      tempo_offset_seconds * out_audio->header().n_samples_per_sec;

  for (int t = 0; t < n_ticks; ++t)
  {
    const int sample_offset =
        tempo_offset_n_samples + (int)(t * tempo_n_samples_per_beat);
    out_audio->OverdubSamples(sample_offset, metronome_tick,
                              1.0 - tick_amplitude);
  }
}

void MakeLevelOneFftProcessorConfig (
    double level_1_n_ffts_per_sec_target,
    int audio_n_samples_per_sec,
    int level_2_fft_spacing,
    int n_total_samples,
    pipeline::FftProcessorConfig* cfg_out)
{
  // Level-2 FFT spacing MUST be even to avoid a future SIMD misalignment when
  // the FftProcessor distributes the real-to-complex transforms among multiple
  // pipeline processing threads.
  // Note that 2*sizeof(double) == 16, and 16 is the SIMD alignment boundary size.
  assert(0 == level_2_fft_spacing % 2);

  // The actual # FFTs per second may differ from the desired value
  int level_1_fft_len =
      static_cast<int>(audio_n_samples_per_sec / level_1_n_ffts_per_sec_target);

  // Shift level_1_fft_len to the nearest value that would yield an output
  // of length divisible by level_2_fft_spacing; note that the type of
  // FFT transform we're doing -- real-to-complex -- outputs fewer values than
  // were input. Have a look at the documentation for
  // FftProcessorConfig::out_len_per_fft() for an explanation of the weird
  // footwork here.
  if (0 != (level_1_fft_len/2+1) % level_2_fft_spacing) {
    level_1_fft_len -= 2 * ((level_1_fft_len/2+1) % level_2_fft_spacing);
  }
  const int min_level_1_fft_len = 2 * (level_2_fft_spacing-1);
  if (level_1_fft_len < min_level_1_fft_len) {
    level_1_fft_len = min_level_1_fft_len;
  }

  const int n_ffts = n_total_samples / level_1_fft_len;

  vector<double> fft_window;
  MakeHanningWindow(level_1_fft_len, &fft_window);

  cfg_out->set_type      (FftType::kRealToComplex);
  cfg_out->set_fft_len   (level_1_fft_len);
  cfg_out->set_n_ffts    (n_ffts);
  cfg_out->set_in_stride (1);
  cfg_out->set_in_dist   (level_1_fft_len);
  cfg_out->set_out_stride(1);
  cfg_out->set_window    (fft_window);

  // The output of a r2c transform is symmetric, and FFTW omits the mirror
  // image:  http://www.fftw.org/fftw3_doc/Real_002ddata-DFT-Array-Format.html
  cfg_out->set_out_dist(level_1_fft_len/2+1);

  cfg_out->set_fftw_planner_time_limit_sec(kDefaultFftwPlannerTimeLimitSec);
}

/** Construct the FftProcessor configuration for the second-level FFTs.
 * The second-level FFT is performed on the output of the first-level FFT.
 * @param level_2_fft_len FFT length.
 * @param level_2_fft_spacing Distance, in bins (first-level FFT output),
 * between consecutive second-level FFTs.
 * @param level_1_cfg The config used for the first level FFTs.
 * @param
 */
void MakeLevelTwoFftProcessorConfig (
    int level_2_fft_len,
    int level_2_fft_spacing,
    const pipeline::FftProcessorConfig& level_1_cfg,
    pipeline::FftProcessorConfig* cfg_out)
{
  // This condition should have been validated prior to pipeline construction
  assert(0 == level_1_cfg.out_len_per_fft() % level_2_fft_spacing);

  const int n_ffts_per_level_1_fft_group =
      level_1_cfg.out_len_per_fft() / level_2_fft_spacing;
  const int level_2_n_ffts =
      n_ffts_per_level_1_fft_group * (level_1_cfg.n_ffts() / level_2_fft_len);

  vector<double> fft_window;
  // Uniform window to disable windowing. If we apply a non-uniform window to
  // the level-2 cornerturn FFT input, we'll likely have at least a slightly
  // adverse impact on the tempo calculation.
  MakeUniformWindow(level_2_fft_len, &fft_window);

  cfg_out->set_type      (FftType::kRealToComplex);
  cfg_out->set_fft_len   (level_2_fft_len);
  cfg_out->set_n_ffts    (level_2_n_ffts);
  cfg_out->set_in_stride (level_1_cfg.out_len_per_fft());
  cfg_out->set_in_dist   (level_2_fft_spacing);
  cfg_out->set_out_stride(1);
  cfg_out->set_window    (fft_window);

  // Important to ensure that FFTW doesn't wreck the input level-1 spectrum
  // power data while calculating the level-2 output data, as we'll need the
  // level-1 powers again later when estimating the first-beat time offset
  cfg_out->set_fftw_planner_flags(FFTW_PRESERVE_INPUT);

  // The output of a r2c transform is symmetric, and FFTW omits the mirror
  // image:  http://www.fftw.org/fftw3_doc/Real_002ddata-DFT-Array-Format.html
  cfg_out->set_out_dist(level_2_fft_len/2+1);
}


int main (int argc, char** argv)
{
  const string progname = argv[0];

  double level_1_n_ffts_per_sec_target = kDefaultLevelOneFftsPerSecond;
  int level_2_fft_spacing = kDefaultLevelTwoFftSpacing;
  double level_2_fft_span_sec = kDefaultLevelTwoFftSpanSeconds;
  int n_pipeline_proc_threads = kDefaultNPipelineProcThreads;

  string in_file_path;

  // Handle the command line arguments.
  char ch = 0;
  while ( EOF != (ch = getopt(argc, argv, "f:hi:t:")) ) {
    switch (ch) {
    case 'f': level_1_n_ffts_per_sec_target = strtof(optarg, nullptr); break;
    case 'h': Usage(progname); return EXIT_SUCCESS;
    case 'i': in_file_path = optarg; break;
    case 't': n_pipeline_proc_threads = strtol(optarg, nullptr, 10); break;
    //case 'l': level_2_fft_len = strtol(optarg, nullptr, 10); break;
    //case 's': level_2_fft_spacing = strtol(optarg, nullptr, 10); break;
    case '?': Usage(progname); return EXIT_FAILURE;
    default:  Usage(progname); return EXIT_FAILURE;
    }
  }

  if (in_file_path.empty()) {
    cerr << "Input .wav file undefined (-i option).\n";
    Usage(progname);
    return EXIT_FAILURE;
  }

  string result;

  cout << "Processing \"" << in_file_path << "\" ...\n";

  unique_ptr<PcmMonoAudioData> in_audio (new PcmMonoAudioData(in_file_path));
  if (!in_audio->init_error().empty()) {
    cerr << "PCM audio file error for '" << in_file_path << "': "
         << in_audio->init_error() << "\n";
    return EXIT_FAILURE;
  }

  const int n_samples_per_sec = in_audio->header().n_samples_per_sec;

  const int n_total_in_samples = in_audio->n_samples();

  auto thread_squad = make_shared<ThreadSquad>(n_pipeline_proc_threads);

  // Note that in_audio ceases to point to a PcmMonoAudioData upon completion
  // of the AudioFileReader construction
  auto audio_file_reader = make_shared<pipeline::AudioFileReader>(
      thread_squad, std::move(in_audio));

  auto type_converter =
      make_shared<pipeline::TypeConverter<int16_t, double>>();

  pipeline::FftProcessorConfig level_1_fft_proc_cfg;
  MakeLevelOneFftProcessorConfig(level_1_n_ffts_per_sec_target,
                                 n_samples_per_sec,
                                 level_2_fft_spacing,
                                 n_total_in_samples,
                                 &level_1_fft_proc_cfg);

  auto level_1_fft_proc =
      make_shared<pipeline::FftProcessor>(level_1_fft_proc_cfg);

  const double level_1_n_ffts_per_sec_effective =
      static_cast<double>(n_samples_per_sec) / level_1_fft_proc_cfg.fft_len();

  const int level_2_fft_len =
      static_cast<int>(level_2_fft_span_sec * level_1_n_ffts_per_sec_effective);

  auto complex_power = make_shared<pipeline::ComplexPower>();

  pipeline::FftProcessorConfig level_2_fft_proc_cfg;
  MakeLevelTwoFftProcessorConfig(level_2_fft_len,
                                 level_2_fft_spacing,
                                 level_1_fft_proc_cfg,
                                 &level_2_fft_proc_cfg);

  auto level_2_fft_proc =
      make_shared<pipeline::FftProcessor>(level_2_fft_proc_cfg);

  // NOTE: This is not a min/max mismatch. Samples per beat is inversely
  // proportional to beats per minute.
  const double min_acceptable_tempo_samples_per_beat =
      BpmToSamplesPerBeat(kMaxAcceptableTempoBpm, n_samples_per_sec);

  // NOTE: This is not a min/max mismatch. Samples per beat is inversely
  // proportional to beats per minute.
  const double max_acceptable_tempo_samples_per_beat =
      BpmToSamplesPerBeat(kMinAcceptableTempoBpm, n_samples_per_sec);

  auto tempo_estimator = make_shared<pipeline::TempoEstimator>(
      level_1_fft_proc_cfg,
      level_2_fft_proc_cfg,
      complex_power,
      min_acceptable_tempo_samples_per_beat,
      max_acceptable_tempo_samples_per_beat );

  // Configure pipeline connections
  audio_file_reader->ConnectOutput(type_converter);
  type_converter   ->ConnectOutput(level_1_fft_proc);
  level_1_fft_proc ->ConnectOutput(complex_power);
  complex_power    ->ConnectOutput(level_2_fft_proc);
  level_2_fft_proc ->ConnectOutput(tempo_estimator);

  sighandler_data_pipeline_head =
      static_pointer_cast<pipeline::AbstractPipelineHead>(audio_file_reader);
  signal(SIGTERM, HandleSignal);
  signal(SIGINT,  HandleSignal);

  shared_ptr<pipeline::AbstractPipelineHead> pipeline_head = audio_file_reader;

  pipeline_head->Start();
  pipeline_head->WaitForProcessingCompletion();

  if (!tempo_estimator->process_error_str().empty()) {
    cerr << "TempoEstimator pipeline node failed: "
         << tempo_estimator->process_error_str() << "\n";
    return EXIT_FAILURE;
  }

  if (!pipeline_head->pipeline_error_str().empty()) {
    cerr << "Pipeline processing failed: "
         << pipeline_head->pipeline_error_str() << "\n";
    return EXIT_FAILURE;
  }

  const double tempo_offset_seconds =
      tempo_estimator->tempo_offset_samples()
      / static_cast<double>(n_samples_per_sec);

  const double average_tempo_bpm =
      n_samples_per_sec * 60.0
      / tempo_estimator->average_tempo_samples_per_beat();

  cout << "Tempo estimate (beats per minute): "
       << average_tempo_bpm
       << ", offset = " << tempo_offset_seconds << " sec"
       << "\n";

  // Now write a new PCM audio file with a metronome sound overdubbed per the
  // tempo beat rate and offset

  shared_ptr<PcmMonoAudioData> out_audio = audio_file_reader->ClonePcmAudioData();

  OverdubMetronome(
      tempo_offset_seconds,
      average_tempo_bpm,
      kMetronomeTickAmplitude,
      kMetronomeToneFreq,
      kMetronomeTickDurationSec,
      out_audio.get() );

  size_t final_slash_index = in_file_path.find_last_of('/');
  const string in_file_name =
      ( string::npos == final_slash_index ?
            in_file_path : in_file_path.substr(final_slash_index+1) );

  const string out_file_path = "./TEMPO_OUTPUT_" + in_file_name;

  cout << "Writing output file \"" << out_file_path << "\" ...\n";

  result = out_audio->WriteToFile(out_file_path);
  if (!result.empty()) {
    cerr << "Failed to write output file with overdubbed metronome '"
         << out_file_path << "': " << result << "\n";
    return EXIT_FAILURE;
  }

  cout << "Done!\n";

  return EXIT_SUCCESS;
}
