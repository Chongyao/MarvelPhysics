file(GLOB header *.h)
file(GLOB source *.cc *.c)
if (${APPLE})
  add_library(FEM SHARED ${source})
else ()
  cuda_add_library(FEM SHARED ${source})
endif (${APPLE})

target_link_libraries(FEM
  marvel-core
 )



