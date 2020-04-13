cmake_minimum_required(VERSION 2.6)
project(mprgp_solver)

# Detects whether this is a top-level project
get_directory_property(MPRGP_SOLVER_PARENT_DIR PARENT_DIRECTORY)
if(NOT MPRGP_SOLVER_PARENT_DIR)
  set(MPRGP_SOLVER_TOPLEVEL_PROJECT ON)
else(NOT MPRGP_SOLVER_PARENT_DIR)
  set(MPRGP_SOLVER_TOPLEVEL_PROJECT OFF)
endif(NOT MPRGP_SOLVER_PARENT_DIR)

# Build tests and samples
option(MPRGP_SOLVER_BUILD_TESTS   "Build mprgp solver test"    ${MPRGP_SOLVER_TOPLEVEL_PROJECT})
option(MPRGP_SOLVER_BUILD_SAMPLES "Build mprgp solver samples" ${MPRGP_SOLVER_TOPLEVEL_PROJECT})

find_package(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
  message(STATUS "EIGEN3_INCLUDE_DIRS: ${EIGEN3_INCLUDE_DIRS}")
  message(STATUS "EIGEN3_VERSION     : ${EIGEN3_VERSION_STRING}")
  include_directories(${EIGEN3_INCLUDE_DIRS})
endif (EIGEN3_FOUND)


# include some directories for better coding
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/utils)

# include the subdirectories
add_subdirectory(src)

# According to the flags, build the test and sampels or not.
# if(MPRGP_SOLVER_BUILD_TESTS)
  add_subdirectory(tests)
# endif(MPRGP_SOLVER_BUILD_TESTS)

if(MPRGP_SOLVER_BUILD_SAMPLES)
  add_subdirectory(samples)
endif(MPRGP_SOLVER_BUILD_SAMPLES)
