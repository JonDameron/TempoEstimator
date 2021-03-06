cmake_minimum_required (VERSION 2.8)

if (NOT DEFINED PROJECT_BASE_DIR)
  message (FATAL_ERROR "PROJECT_BASE_DIR not defined; set it in the toplevel CMakeLists")
endif (NOT DEFINED PROJECT_BASE_DIR)

if (NOT DEFINED RUNTIME_INSTALL_PATH)
  message (FATAL_ERROR "RUNTIME_INSTALL_PATH not defined; set it in the toplevel CMakeLists")
endif (NOT DEFINED RUNTIME_INSTALL_PATH)

if (NOT DEFINED INSTALL_PERMISSIONS)
  message (FATAL_ERROR "INSTALL_PERMISSIONS not defined; set it in the toplevel CMakeLists")
endif (NOT DEFINED INSTALL_PERMISSIONS)

if (NOT DEFINED MAX_ERRORS)
  set(MAX_ERRORS 5)
endif(NOT DEFINED MAX_ERRORS)

set (SOURCES
  abstract_pipeline_head.cpp
  abstract_pipeline_node.cpp
  audio_file_reader.cpp
  complex_power.cpp
  fft_processor.cpp
  fft_processor_config.cpp
  tempo_estimator.cpp
  )

add_library (pipeline  STATIC ${SOURCES})

target_link_libraries (pipeline  boost_system boost_thread fftw3 m)
