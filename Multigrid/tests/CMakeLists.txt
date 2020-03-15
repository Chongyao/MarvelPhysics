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

add_executable(V_cycle V_cycle.cc)
target_link_libraries(V_cycle
  marvel-core
  MULTIGRID
  FEM
  )

add_executable(test_GS test_GS.cc)
target_link_libraries(test_GS
  marvel-core
  MULTIGRID
  FEM
  )

add_executable(test_permutation test_permutation.cc)
target_link_libraries(test_permutation
  marvel-core
  MULTIGRID
  FEM
  )


