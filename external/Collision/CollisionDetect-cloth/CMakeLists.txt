file(GLOB coll_cloth_header src/*.h)
file(GLOB coll_cloth_source src/*.cc src/*.c src/*.cpp)
file(GLOB coll_cloth_cuda src/*.cuh src/*.cu)

#find_package(CUDA)
FIND_PACKAGE(CUDA ${CUDA_VERSION} REQUIRED)
STRING (TOLOWER ${CMAKE_SYSTEM_NAME} LOWER_SYSTEM_NAME)
SET(CUDA_TARGET_INCLUDE  ${CUDA_TOOLKIT_ROOT_DIR}-${CUDA_VERSION}/targets/${CMAKE_HOST_SYSTEM_PROCESSOR}-${LOWER_SYSTEM_NAME}/include)

#TODO: make this available for others
SET(CUDA_NVCC_FLAGS --compiler-options "-std=c++11")
SET(CUDA_HELPER_PATH ${PROJECT_SOURCE_DIR}/external/cuda-samples/Common)

include_directories(${CUDA_HELPER_PATH})

CUDA_ADD_LIBRARY(coll_cloth ${coll_cloth_header} ${coll_cloth_source} ${coll_cloth_cuda})
target_link_libraries(coll_cloth
  GLU glut GL cublas
  ${CUDA_cusparse_LIBRARY}
)
