cmake_minimum_required (VERSION 2.8)

set (SOURCES
  main.cpp
  )

add_executable (tempo_estimator  ${SOURCES})

target_link_libraries (tempo_estimator  pipeline m boost_thread fftw3)
