cmake_minimum_required(VERSION 2.6)
project(chaos_utils)

# Detects whether this is a top-level project
get_directory_property(CHAOS_UTILS_PARENT_DIR PARENT_DIRECTORY)
if(NOT CHAOS_UTILS_PARENT_DIR)
  set(CHAOS_UTILS_TOPLEVEL_PROJECT ON)
else(NOT CHAOS_UTILS_PARENT_DIR)
  set(CHAOS_UTILS_TOPLEVEL_PROJECT OFF)
endif(NOT CHAOS_UTILS_PARENT_DIR)

# Build tests and samples
option(CHAOS_UTILS_BUILD_TESTS   "Build chaos utils test"    ${CHAOS_UTILS_TOPLEVEL_PROJECT})
option(CHAOS_UTILS_BUILD_SAMPLES "Build chaos utils samples" ${CHAOS_UTILS_TOPLEVEL_PROJECT})

# CHAOS_UTILS_LIBRARY_FLAG speeds up the generation of multiple binaries
# at the cost of a longer initial compilation time
# (by default, chaos_utils_library is off since we can use this utils as
# header only library.)
option(CHAOS_UTILS_LIBRARY_FLAG "Use chaos utils as library" OFF)

# Some components of chaos utils library flags
# (by default, they are all off since we can use them as header only library)
option(CHAOS_LOG_LIBRARY_FLAG   "Use logger util as library" OFF)
option(CHAOS_TIMER_LIBRARY_FLAG "Use timer  util as library" OFF)

# add some compile definition to the source files according to the options
# the global chaos utils library flags
if(CHAOS_UTILS_LIBRARY_FLAG)
  add_definitions(-DCHAOS_UTILS_LIBRARY_FLAG)
endif(CHAOS_UTILS_LIBRARY_FLAG)

# the logger library flag
if(CHAOS_LOG_LIBRARY_FLAG)
  add_definitions(-DCHAOS_LOG_LIBRARY_FLAG)
endif(CHAOS_LOG_LIBRARY_FLAG)

# the timer library flag
if(CHAOS_TIMER_LIBRARY_FLAG)
  add_definitions(-DCHAOS_TIMER_LIBRARY_FLAG)
endif(CHAOS_TIMER_LIBRARY_FLAG)

# include some directories for better coding
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/logger)
include_directories(${PROJECT_SOURCE_DIR}/timer)
include_directories(${PROJECT_SOURCE_DIR}/oshelper)

# include the subdirectories
add_subdirectory(logger)
add_subdirectory(timer)

# According to the flags, build the test and sampels or not.
if(CHAOS_UTILS_BUILD_TESTS)
  add_subdirectory(tests)
endif(CHAOS_UTILS_BUILD_TESTS)

if(CHAOS_UTILS_BUILD_SAMPLES)
  add_subdirectory(samples)
endif(CHAOS_UTILS_BUILD_SAMPLES)
