add_executable(test_elas test_elas.cc)
target_link_libraries(test_elas
  marvel-core
  FEM)

add_executable(test_three_component test_three_component.cc)
target_link_libraries(test_three_component
  marvel-core
  FEM)
