project(MarvalPhysics)
cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_CXX_FLAGS "-fopenmp -std=c++11 -fpermissive")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${CMAKE_MODULE_PATH}")

include(geo_sim_sdk)
include_geo_sim_sdk()
link_geo_sim_sdk()
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

# PETSC
set(PETSC_DIR "/usr/lib/petscdir/3.7.7")
set(PETSC_ARCH "x86_64-linux-gnu-real")
find_package(PETSc)
if(PETSC_FOUND)
  include_directories(/usr/lib/openmpi/include)
  include_directories(/usr/lib/openmpi/lib)
  message("-- PETSc @ ${PETSC_INCLUDE_DIR}")
  message("-- PETSc version: ${PETSC_VERSION}")
  include_directories(${PETSC_INCLUDE_DIR})
endif()

# UMFPACK
find_package(UMFPACK REQUIRED)
if(UMFPACK_FOUND)
  message("-- Umfpack @ ${UMFPACK_INCLUDES}")
  include_directories(${UMFPACK_INCLUDES})
endif(UMFPACK_FOUND)

# SUITESPARSE
find_package(SuiteSparse REQUIRED)
if(SuiteSparse_FOUND)
  message("-- SuiteSparse @ ${SuiteSparse_INCLUDE_DIRS}")
  message("-- SuiteSparse libs @ ${SuiteSparse_LIBRARIES}")
  include_directories(${SuiteSparse_INCLUDE_DIRS})
endif(SuiteSparse_FOUND)

# LBFGS
find_package(LBFGS REQUIRED)
if(LBFGS_FOUND)
  message("-- LBFGS @ ${LBFGS_INCLUDE_DIRS}")
  include_directories(${LBFGS_INCLUDE_DIRS})
endif(LBFGS_FOUND)

# IPOPT
find_package(IPOPT REQUIRED)
if(IPOPT_FOUND)
  message("-- IPOPT @ ${IPOPT_INCLUDE_DIRS}")
  message("-- IPOPT lib @ ${IPOPT_LIBRARIES}")
  include_directories(${IPOPT_INCLUDE_DIRS})
endif(IPOPT_FOUND)

# ALGLIB
find_package(ALGLIB REQUIRED)
if(ALGLIB_FOUND)
  message("-- Alglib @ ${ALGLIB_INCLUDES}")
  message("-- Alglibs @ ${ALGLIB_LIBRARIES}")
  include_directories(${ALGLIB_INCLUDES})
endif(ALGLIB_FOUND)

# ARPACK
find_package(Arpack REQUIRED)

# GLOG
find_package(Glog REQUIRED)
if(GLOG_FOUND)
  message("-- GLOG @ ${GLOG_INCLUDE_DIRS}")
  include_directories(${GLOG_INCLUDE_DIRS})
endif(GLOG_FOUND)

# GTEST
find_package(GTest REQUIRED)
if(GTEST_FOUND)
  message("-- GTEST @ ${GTEST_INCLUDE_DIRS}")
  include_directories(${GTEST_INCLUDE_DIRS})
endif(GTEST_FOUND)

# JSONCPP
find_package(JsonCpp REQUIRED)
if(JsonCpp_FOUND)
  message("-- JsonCpp @ ${JsonCpp_INCLUDE_DIRS}")
  include_directories(${JsonCpp_INCLUDE_DIRS})
endif(JsonCpp_FOUND)

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




add_subdirectory(external)
# add_subdirectory(examples)
# add_subdirectory(test)
# add_subdirectory(coarsen)
add_subdirectory(Point_Sys)
add_subdirectory(src)

