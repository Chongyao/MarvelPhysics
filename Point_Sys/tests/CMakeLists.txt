add_executable(test_genpoints test_genpoints.cc)
target_link_libraries(test_genpoints
  bigbang
  PointSys
  zjucad-ptree
  jtf-mesh)

add_executable(test_wdnb test_wdnb.cc)
target_link_libraries(test_wdnb
  bigbang
  # PointSys
  tbb
  zjucad-ptree
  jtf-mesh)
