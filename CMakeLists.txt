cmake_minimum_required (VERSION 2.8)

project (tempo-estimator)

if (NOT DEFINED EXTRA_C_DEFINES)
  set (EXTRA_C_DEFINES "")
endif (NOT DEFINED EXTRA_C_DEFINES)

set (CMAKE_VERBOSE_MAKEFILE ON)

set (PROJECT_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

set (CMAKE_BUILD_TYPE Debug)

set (ADDITIONAL_CXX_FLAGS "${EXTRA_C_DEFINES}")

set (LIBRARY_INSTALL_PATH "${PROJECT_BASE_DIR}/lib")
set (RUNTIME_INSTALL_PATH "${PROJECT_BASE_DIR}/bin")

set (INSTALL_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Where to put the newly created binaries and libraries
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

set (CMAKE_INSTALL_RPATH ${LIBRARY_INSTALL_PATH})
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
set (CMAKE_SKIP_BUILD_RPATH false)

# Using fmax-errors option to avoid excessive compiler output
set (CMAKE_CXX_FLAGS
  "-std=c++0x -fmax-errors=5 -Wall ${ADDITIONAL_CXX_FLAGS} ")

set (CMAKE_CXX_FLAGS_DEBUG
  "-ggdb3")

set (CMAKE_CXX_FLAGS_RELEASE
  "-O3")

include_directories ("${PROJECT_BASE_DIR}")

add_subdirectory (pipeline)
add_subdirectory (tempo-estimator)
add_subdirectory (utility)

target_link_libraries (tempo_estimator  pipeline m boost_thread fftw3)

SET(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
SET(CPACK_DEBIAN_PACKAGE_SECTION "Science")
SET(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSING})
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "jon.g.dameron@gmail.com")
SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION  "Tempo estimation for musical audio data")
SET(CPACK_PACKAGE_VERSION "1-0-0")
SET(CPACK_PACKAGE_FILE_NAME "tempo_estimator_${CPACK_PACKAGE_VERSION}")
SET(CPACK_GENERATOR "DEB")
# SET(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_CURRENT_SOURCE_DIR}/release/postinst")

include(CPack)
