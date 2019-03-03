file(GLOB header *.h)
file(GLOB source *.cc *.c)
cuda_add_library(PointSys SHARED ${source})
target_link_libraries(PointSys
  coll_cloth
  # WindingNumber
  tbb
#  marvel-core
  )
