cmake_minimum_required(VERSION 2.6)

file(GLOB coll_cloth_header src/*.h)
file(GLOB coll_cloth_source src/*.cc src/*.c src/*.cpp)
file(GLOB coll_cloth_cuda src/*.cuh src/*.cu)

#find_package(CUDA)
FIND_PACKAGE(CUDA ${CUDA_VERSION} REQUIRED)
SET(CUDA_NVCC_FLAGS --compiler-options "-std=c++11")
SET(CUDA_HELPER_PATH ${PROJECT_SOURCE_DIR}/external/cuda-samples/Common)

include_directories(${CUDA_HELPER_PATH})

CUDA_ADD_LIBRARY(coll_v3 ${coll_cloth_header} ${coll_cloth_source} ${coll_cloth_cuda})
target_link_libraries(coll_v3
  GLU glut GL
  ${CUDA_cusparse_LIBRARY}
)
