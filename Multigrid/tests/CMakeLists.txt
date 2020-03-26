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
  ${Boost_LIBRARIES}    
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

add_executable(test_sparsification test_sparsification.cc)
target_link_libraries(test_sparsification
  marvel-core
  MULTIGRID
  FEM
  )

add_executable(test_spectra test_spectra.cc)



