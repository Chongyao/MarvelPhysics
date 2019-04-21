project(MarvalPhysics)
cmake_minimum_required(VERSION 2.8)

# list(APPEND CMAKE_CXX_FLAGS "-fopenmp -std=c++11 -O3 -fpermissive -march=native")
list(APPEND CMAKE_CXX_FLAGS "-fopenmp -std=c++11 -O0 -fpermissive -march=native")

# set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${CMAKE_MODULE_PATH}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# include(geo_sim_sdk)
# include_geo_sim_sdk()
# link_geo_sim_sdk()
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/external)

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

#CUDA
FIND_PACKAGE(CUDA ${CUDA_VERSION} REQUIRED)
STRING (TOLOWER ${CMAKE_SYSTEM_NAME} LOWER_SYSTEM_NAME)
SET(CUDA_TARGET_INCLUDE  ${CUDA_TOOLKIT_ROOT_DIR}-${CUDA_VERSION}/targets/${CMAKE_HOST_SYSTEM_PROCESSOR}-${LOWER_SYSTEM_NAME}/include)
message(${CUDA_INCLUDE_DIRS})
message(${CUDA_TOOLKIT_ROOT_DIR})
#TODO: make this available for others
SET(CUDA_CUBLAS_LIBRARIES /usr/local/cuda/lib64/libcublas.so)
SET(CUDA_NVCC_FLAGS --compiler-options "-std=c++11")
SET(CUDA_HELPER_PATH ${CUDA_TOOLKIT_ROOT_DIR}/samples/common/inc)
include_directories(${CUDA_HELPER_PATH})
list(APPEND CUDA_NVCC_FLAGS "-g -G")

# # PETSC
# set(PETSC_DIR "/usr/lib/petscdir/3.7.7")
# set(PETSC_ARCH "x86_64-linux-gnu-real")
# find_package(PETSc)
# if(PETSC_FOUND)
#   include_directories(/usr/lib/openmpi/include)
#   include_directories(/usr/lib/openmpi/lib)
#   message("-- PETSc @ ${PETSC_INCLUDE_DIR}")
#   message("-- PETSc version: ${PETSC_VERSION}")
#   include_directories(${PETSC_INCLUDE_DIR})
# endif()

# UMFPACK
# find_package(UMFPACK REQUIRED)
# if(UMFPACK_FOUND)
#   message("-- Umfpack @ ${UMFPACK_INCLUDES}")
#   include_directories(${UMFPACK_INCLUDES})
# endif(UMFPACK_FOUND)

# SUITESPARSE
# find_package(SuiteSparse REQUIRED)
# if(SuiteSparse_FOUND)
#   message("-- SuiteSparse @ ${SuiteSparse_INCLUDE_DIRS}")
#   message("-- SuiteSparse libs @ ${SuiteSparse_LIBRARIES}")
#   include_directories(${SuiteSparse_INCLUDE_DIRS})
# endif(SuiteSparse_FOUND)

# LBFGS
# find_package(LBFGS REQUIRED)
# if(LBFGS_FOUND)
#   message("-- LBFGS @ ${LBFGS_INCLUDE_DIRS}")
#   include_directories(${LBFGS_INCLUDE_DIRS})
# endif(LBFGS_FOUND)

# IPOPT
# find_package(IPOPT REQUIRED)
# if(IPOPT_FOUND)
#   message("-- IPOPT @ ${IPOPT_INCLUDE_DIRS}")
#   message("-- IPOPT lib @ ${IPOPT_LIBRARIES}")
#   include_directories(${IPOPT_INCLUDE_DIRS})
# endif(IPOPT_FOUND)

# ALGLIB
# find_package(ALGLIB REQUIRED)
# if(ALGLIB_FOUND)
#   message("-- Alglib @ ${ALGLIB_INCLUDES}")
#   message("-- Alglibs @ ${ALGLIB_LIBRARIES}")
#   include_directories(${ALGLIB_INCLUDES})
# endif(ALGLIB_FOUND)

# # ARPACK
# find_package(Arpack REQUIRED)

# # GLOG
# find_package(Glog REQUIRED)
# if(GLOG_FOUND)
#   message("-- GLOG @ ${GLOG_INCLUDE_DIRS}")
#   include_directories(${GLOG_INCLUDE_DIRS})
# endif(GLOG_FOUND)

# # GTEST
# find_package(GTest REQUIRED)
# if(GTEST_FOUND)
#   message("-- GTEST @ ${GTEST_INCLUDE_DIRS}")
#   include_directories(${GTEST_INCLUDE_DIRS})
# endif(GTEST_FOUND)

# # JSONCPP
# find_package(JsonCpp REQUIRED)
# if(JsonCpp_FOUND)
#   message("-- JsonCpp @ ${JsonCpp_INCLUDE_DIRS}")
#   include_directories(${JsonCpp_INCLUDE_DIRS})
# endif(JsonCpp_FOUND)

#TBB
find_package(TBB REQUIRED)
if(TBB_FOUND)
  message("-- TBB @ ${TBB_INCLUDE_DIRS}")
  include_directories(${TBB_INCLUDE_DIRS})
endif(TBB_FOUND)
# CUDA
# find_package(CUDA)		
# if(CUDA_FOUND)
#   add_definitions(-DUSE_CUDA)
#   message("-- CUDA @ ${CUDA_INCLUDE_DIRS}")
#   include_directories(${CUDA_INCLUDE_DIRS})
# endif(CUDA_FOUND)

#IGL
# option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
# option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
# option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
# find_package(LIBIGL REQUIRED QUIET)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# find_package(OpenMP REQUIRED)
# if(NOT TARGET OpenMP::OpenMP_CXX)
#     add_library(OpenMP_TARGET INTERFACE)
#     add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
#     target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
#     find_package(Threads REQUIRED)
#     target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
#     target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
# endif()


add_subdirectory(external)
add_subdirectory(Point_Sys)
add_subdirectory(src)
# add_subdirectory(test)

