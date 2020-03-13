add_executable(test_SP_COMP test_SP_COMP.cc)
target_link_libraries(test_SP_COMP
  marvel-core
  MULTIGRID
  )

add_executable(test_inverse_map test_inverse_map.cc)
target_link_libraries(test_inverse_map
  marvel-core
  MULTIGRID
  FEM
  )


