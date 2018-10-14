add_executable(test_genpoints test_genpoints.cc)
target_link_libraries(test_genpoints
  PointSys
  zjucad-ptree
  jtf-mesh)
