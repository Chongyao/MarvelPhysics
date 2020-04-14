file(GLOB header *.h)
file(GLOB source *.cc *.c)


if (${APPLE})
  add_library(MULTIGRID SHARED ${source})
else ()
  cuda_add_library(MULTIGRID SHARED ${source})
endif (${APPLE})

target_link_libraries(MULTIGRID
  marvel-core
  )
