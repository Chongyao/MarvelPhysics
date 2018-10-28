file(GLOB header *.h)
file(GLOB source *.cc *.c)
add_library(PointSys SHARED ${source})
target_link_libraries(PointSys
  # WindingNumber
  tbb
#  marvel-core
  )
