file(GLOB header *.h)
file(GLOB source *.cc *.c)
cuda_add_library(FEM SHARED ${source})
target_link_libraries(FEM
  marvel-core
 )



