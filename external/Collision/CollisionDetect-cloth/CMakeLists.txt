file(GLOB coll_cloth_header src/*.h)
file(GLOB coll_cloth_source src/*.cc src/*.c src/*.cpp)
file(GLOB coll_cloth_cuda src/*.cuh src/*.cu)


CUDA_ADD_LIBRARY(coll_cloth ${coll_cloth_header} ${coll_cloth_source} ${coll_cloth_cuda})
target_link_libraries(coll_cloth
  GLU glut GL
  ${CUDA_cusparse_LIBRARY}
  ${CUDA_CUBLAS_LIBRARIES}
)
