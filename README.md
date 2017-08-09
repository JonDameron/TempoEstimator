Musical Tempo Estimation Utility and Associated Digital Signal Processing Libraries  
(c) 2017 Jonathan G. Dameron | jon.g.dameron@gmail.com

C++11 utility that performs multi-threaded digital signal processing of raw audio data to estimate the musical tempo.

Build Instructions
------------------

Required packages:

libfftw3-dev  
libboost-system-dev  
libboost-thread-dev

    cd <source directory>
    mkdir build
    cd build
    cmake ..
    make

If the build succeeds, the output binary will appear in "bin".

Usage Instructions
------------------

    bin/tempo_estimator -i INPUT_WAV_FILE [other options]

Run with -h to see the full usage message with a list of all options.

The Audio File Tempo Estimator reads music data in the form of a .wav file (uncompressed PCM), attempts to estimate the tempo of the music, and outputs a new .wav file with two alternating metronome click sounds. One marks the beats, the other (quieter) marks the offbeats.

The utility usually produces reasonably accurate output for music files that have a rigid tempo with only small fluctuations; if the song was originally recorded with a metronome for reference, this condition is probably met. Tracking of significant tempo fluctuations over the course of a song is not currently supported.

The output file name will match that of the input, with extension ".TEMPO.wav" instead of ".wav".

Sample audio input files are available in the "sample_input" directory. All provided songs are in the US public domain. See included text file in "sample_input".

Source Code Notes
-----------------

For sake of simplicity, most functions that can fail return their error status in the form of a string; an empty string is indicative of success.

