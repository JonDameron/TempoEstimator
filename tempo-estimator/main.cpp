#include <csignal>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "utility/pcm_mono_audio_data.hpp"
#include "utility/util.hpp"
#include "pipeline/audio_file_reader.hpp"
#include "pipeline/fft_processor.hpp"
#include "pipeline/fft_processor_config.hpp"
#include "pipeline/tempo_estimator.hpp"
#include "pipeline/type_converter.hpp"

using namespace std;

/** Default number of FFTs to execute per second while processing the raw
 * audio data.  Setting this too low will reduce the accuracy of the tempo
 * calculation; 200 < value < (sampleRate / 50) will probably generate a
 * reasonably accurate estimate.
 */
static const int    kDefaultLevelOneFftsPerSecond = 200;
static const int    kDefaultLevelTwoFftLen = kDefaultLevelOneFftsPerSecond * 30;
static const int    kDefaultLevelTwoFftSpacing = 10;
static const double kDefaultFftwPlannerTimeLimitSec = 1.0;

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
       << "  -f  Frequency profile granularity (N per second) [default = "
           << kDefaultLevelOneFftsPerSecond << "]\n"
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

void MakeLevelOneFftProcessorConfig (
    int level_1_n_ffts_per_sec,
    int audio_n_samples_per_sec,
    int level_2_fft_spacing,
    int n_total_samples,
    pipeline::FftProcessorConfig* cfg_out)
{
  int level_1_fft_len = audio_n_samples_per_sec / level_1_n_ffts_per_sec;
  // Round up to nearest multiple of level_2_fft_spacing
  level_1_fft_len += level_1_fft_len % level_2_fft_spacing;

  const int n_ffts  = n_total_samples / level_1_fft_len;

  vector<double> fft_window;
  MakeHanningWindow(level_1_fft_len, &fft_window);

  cfg_out->set_fft_len   (level_1_fft_len);
  cfg_out->set_n_ffts    (n_ffts);
  cfg_out->set_in_stride (1);
  cfg_out->set_in_dist   (level_1_fft_len);
  cfg_out->set_out_stride(1);
  cfg_out->set_out_dist  (level_1_fft_len);
  cfg_out->set_window    (fft_window);

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
  assert(0 == level_1_cfg.fft_len() % level_2_fft_spacing);

  const int n_ffts_per_level_1_fft_group =
      level_1_cfg.fft_len() / level_2_fft_spacing;
  const int level_2_n_ffts =
      n_ffts_per_level_1_fft_group * (level_1_cfg.n_ffts() / level_2_fft_len);

  vector<double> fft_window;
  MakeHanningWindow(level_2_fft_len, &fft_window);

  cfg_out->set_fft_len   (level_2_fft_len);
  cfg_out->set_n_ffts    (level_2_n_ffts);
  cfg_out->set_in_stride (level_1_cfg.fft_len());
  cfg_out->set_in_dist   (level_2_fft_spacing);
  cfg_out->set_out_stride(1);
  cfg_out->set_out_dist  (level_2_fft_len);
  cfg_out->set_window    (fft_window);
}

int main (int argc, char** argv)
{
  const string progname = argv[0];

  int level_1_n_ffts_per_sec = kDefaultLevelOneFftsPerSecond;
  int level_2_fft_len = kDefaultLevelTwoFftLen;
  int level_2_fft_spacing = kDefaultLevelTwoFftSpacing;
  string input_file_path;

  // Handle the command line arguments.
  char ch = 0;
  while ( EOF != (ch = getopt(argc, argv, "f:h")) ) {
    switch (ch) {
    case 'f': level_1_n_ffts_per_sec = strtol(optarg, NULL, 10); break;
    case 'h': Usage(progname); return EXIT_SUCCESS;
    case 'i': input_file_path = optarg; break;
    //case 'l': level_2_fft_len = strtol(optarg, NULL, 10); break;
    //case 's': level_2_fft_spacing = strtol(optarg, NULL, 10); break;
    case '?': Usage(progname); return EXIT_FAILURE;
    default:  Usage(progname); return EXIT_FAILURE;
    }
  }

  if (input_file_path.empty()) {
    cerr << "Input .wav file undefined (-i option).\n";
    Usage(progname);
    return EXIT_FAILURE;
  }

  string result;

  // Avoid future SIMD-alignment problems in FftProcessor by simply ensuring
  // that that the total # PCM samples is even.
  const int required_n_samples_divisor = 2;

  PcmMonoAudioData audio_data (input_file_path, required_n_samples_divisor);
  if (!audio_data.init_error().empty()) {
    cerr << "PCM audio file error for '" << input_file_path << "': "
         << audio_data.init_error() << "\n";
    return EXIT_FAILURE;
  }

  const int n_samples_per_sec = audio_data.header().n_samples_per_sec;
  const int n_total_samples   = audio_data.n_samples();

  auto audio_file_reader =
      make_shared<pipeline::AudioFileReader>(audio_data);

  auto type_converter =
      make_shared<pipeline::TypeConverter<int16_t, double>>();

  pipeline::FftProcessorConfig level_1_fft_processor_cfg;
  MakeLevelOneFftProcessorConfig(level_1_n_ffts_per_sec,
                                 audio_data.header().n_samples_per_sec,
                                 level_2_fft_spacing,
                                 audio_data.n_samples(),
                                 &level_1_fft_processor_cfg);

  auto level_1_fft_processor =
      make_shared<pipeline::FftProcessor>(level_1_fft_processor_cfg);

  pipeline::FftProcessorConfig level_2_fft_processor_cfg;
  MakeLevelTwoFftProcessorConfig(level_2_fft_len,
                                 level_2_fft_spacing,
                                 level_1_fft_processor_cfg,
                                 &level_2_fft_processor_cfg);

  auto level_2_fft_processor =
      make_shared<pipeline::FftProcessor>(level_2_fft_processor_cfg);

  //auto tempo_estimator = make_shared<pipeline::TempoEstimator>();

  audio_file_reader->ConnectOutput(type_converter);

  type_converter->ConnectOutput(level_1_fft_processor);

  level_1_fft_processor->ConnectOutput(level_2_fft_processor);

  //level_2_fft_processor->ConnectOutput(tempo_estimator);

  sighandler_data_pipeline_head =
      static_pointer_cast<pipeline::AbstractPipelineHead>(audio_file_reader);
  signal(SIGTERM, HandleSignal);
  signal(SIGINT,  HandleSignal);

  shared_ptr<pipeline::AbstractPipelineHead> pipeline_head = audio_file_reader;

  pipeline_head->Start();
  //tempo_estimator->WaitForEndOfInput();

  PcmMonoAudioData

  cout << "Tempo estimate (beats per minute): "
       << tempo_estimator->average_tempo_bpm() << "\n";

  return EXIT_SUCCESS;
}
