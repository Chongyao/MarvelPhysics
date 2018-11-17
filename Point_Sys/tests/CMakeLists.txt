# add_executable(test_genpoints test_genpoints.cc)
# target_link_libraries(test_genpoints

#   PointSys
#   zjucad-ptree
#   )

# add_executable(test_wdnb test_wdnb.cc)
# target_link_libraries(test_wdnb
#   bigbang
#   # PointSys
#   tbb
#   zjucad-ptree
#   jtf-mesh)


# add_executable(test_energy test_energy.cc)
# target_link_libraries(test_energy
#   PointSys
#   marvel-core
#   ${Boost_LIBRARIES}
#   )

# add_executable(test_linearInterp test_linearInterp.cc)
# target_link_libraries(test_linearInterp
#   PointSys
#   marvel-core
#   ${Boost_LIBRARIES}
#   )


add_executable(test_tet test_tet.cc)
target_link_libraries(test_tet
  PointSys
  marvel-core
  ${Boost_LIBRARIES}
  )

add_executable(test_autodiff test_autodiff.cc)
