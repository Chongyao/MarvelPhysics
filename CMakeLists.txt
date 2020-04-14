project(MarvalPhysics)
cmake_minimum_required(VERSION 2.8)

# list(APPEND CMAKE_CXX_FLAGS "-fopenmp -std=c++11 -O3 -fpermissive -march=native")
if (${APPLE})
  list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -fpermissive -march=native -ftemplate-backtrace-limit=0 -fPIC")
else ()
  list(APPEND CMAKE_CXX_FLAGS "-fopenmp -std=c++11 -O3 -fpermissive -march=native -ftemplate-backtrace-limit=0 -fPIC")
endif (${APPLE})

# set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${CMAKE_MODULE_PATH}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/external)
include_directories(${PROJECT_SOURCE_DIR}/mprgp_solver/src/utils)

# BOOST
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost COMPONENTS filesystem system program_options REQUIRED)
if(Boost_FOUND)
  message("-- Boost @ ${Boost_INCLUDE_DIRS}")
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)
# EIGEN3
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  message("-- Eigen @ ${EIGEN3_INCLUDE_DIR}")
  message("-- Eigen version: ${EIGEN3_VERSION}")
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

if (${APPLE})
  message ("MAC OS X, avoid adding cuda")
else ()
  # CUDA
  FIND_PACKAGE(CUDA ${CUDA_VERSION} REQUIRED)
  STRING (TOLOWER ${CMAKE_SYSTEM_NAME} LOWER_SYSTEM_NAME)
  SET(CUDA_TARGET_INCLUDE  ${CUDA_TOOLKIT_ROOT_DIR}-${CUDA_VERSION}/targets/${CMAKE_HOST_SYSTEM_PROCESSOR}-${LOWER_SYSTEM_NAME}/include)
  #TODO: make this available for others
  SET(CUDA_CUBLAS_LIBRARIES /usr/local/cuda/lib64/libcublas.so)
  SET(CUDA_NVCC_FLAGS --compiler-options "-std=c++11")
  SET(CUDA_HELPER_PATH ${CUDA_TOOLKIT_ROOT_DIR}/samples/common/inc)
  include_directories(${CUDA_HELPER_PATH})
  list(APPEND CUDA_NVCC_FLAGS "-g -G")
endif (${APPLE})

#TBB
find_package(TBB REQUIRED)
if(TBB_FOUND)
  message("-- TBB @ ${TBB_INCLUDE_DIRS}")

  include_directories(${TBB_INCLUDE_DIRS})
endif(TBB_FOUND)

if (NOT ${APPLE})
  find_package(OpenMP REQUIRED)
  if(NOT TARGET OpenMP::OpenMP_CXX)
      add_library(OpenMP_TARGET INTERFACE)
      add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
      target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
      find_package(Threads REQUIRED)
      target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
      target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})

      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif (NOT ${APPLE})

if (${APPLE})
  add_subdirectory(FEM)
  add_subdirectory(src)
  add_subdirectory(Multigrid)


else ()
  add_subdirectory(external)
  # add_subdirectory(Point_Sys)
  add_subdirectory(FEM)
  add_subdirectory(src)
  add_subdirectory(mprgp_solver)
  add_subdirectory(App)


endif (${APPLE})
