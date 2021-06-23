# add_executable(test_data_str test_data_str.cc)
# target_link_libraries(test_data_str
#   marvel-core
#   coll_rigid
#   )
# add_executable(test_openmp test_openmp.cc)


# add_executable(test_bay test_bay.cc)
# cuda_add_executable(test_coll_1 test_coll1.cc)
# target_link_libraries(test_coll_1
#   coll_cloth
#   marvel-core
#   PointSys
#   ${Boost_LIBRARIES}
#   )

# cuda_add_executable(test_coll test_coll.cc)
# target_link_libraries(test_coll
#   coll_rigid
#   marvel-core
#   PointSys
#   ${Boost_LIBRARIES}
#   )
# cuda_add_executable(test_coll_rky test_coll_rky.cc)
# target_link_libraries(test_coll_rky
#   coll_cloth
#   marvel-core
#   PointSys
#   ${Boost_LIBRARIES}
#   )



add_executable(test_genpoints test_genpoints.cc)
target_link_libraries(test_genpoints
  marvel-core
  PointSys
  ${Boost_LIBRARIES}
  )

# add_executable(test_wdnb test_wdnb.cc)
# target_link_libraries(test_wdnb
#   bigbang
#   # PointSys
#   tbb
#   zjucad-ptree
#   jtf-mesh)


# add_executable(test_energy_implicit test_energy_implicit.cc)
# target_link_libraries(test_energy_implicit
#   PointSys
#   marvel-core
#   ${Boost_LIBRARIES}
#   )

# add_executable(test_new_struct test_new_struct.cc)
# target_link_libraries(test_new_struct
#   PointSys
#   marvel-core
#   ${Boost_LIBRARIES}
#   )

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

#add_executable(test_autodiff test_autodiff.cc)
# add_executable(test_SparseMatrix test_SparseMatrix.cc)

# add_executable(blender_meshless blender_meshless.cpp)
# target_link_libraries(blender_meshless
#   coll_cloth
#   PointSys
#   marvel-core
#   ${Boost_LIBRARIES}
#   )

# cuda_add_executable(blender_meshless_coll blender_coll.cpp)
# target_link_libraries(blender_meshless_coll
#   marvel-core
#   coll_rigid
#   PointSys

#   ${Boost_LIBRARIES}
#   )

