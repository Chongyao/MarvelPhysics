
add_executable(test_vtk2surf test_vtk2surf.cc)
target_link_libraries(test_vtk2surf
  marvel-core
  ${Boost_LIBRARIES}
  )

