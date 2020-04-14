SET(CMAKE_CXX_FLAGS "-std=c++0x -fopenmp")

add_executable(example_based example_based.cc)
target_link_libraries(example_based
  marvel-core
  FEM)

add_executable(gen_pose gen_pose.cc)
target_link_libraries(gen_pose
  marvel-core
  FEM)

add_executable(test_elas test_elas.cc)
target_link_libraries(test_elas
  marvel-core
  FEM
  blas
  lapack)

add_executable(test_poisson test_poisson.cc)
target_link_libraries(test_poisson
  marvel-core
  FEM
  blas
  lapack)


SET(CMAKE_CXX_FLAGS "-std=c++0x")


add_executable(test_interp test_interp.cc)
target_link_libraries(test_interp
  marvel-core
  FEM
  blas
  lapack)


add_executable(test_three_component test_three_component.cc)
target_link_libraries(test_three_component
  marvel-core
  FEM)

add_executable(test_extract_surface test_extract_surface.cc)
target_link_libraries(test_extract_surface
  marvel-core
  FEM
  blas
  lapack)
