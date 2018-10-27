enable_language(Fortran)

file(GLOB header *.h)
file(GLOB source *.cc *.c)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
message("this is cmkae current source dir ${CMAKE_CURRENT_SOURCE_DIR}  ")
add_custom_command(OUTPUT potential.f90
    COMMAND "maxima" ARGS  --userdir=${SRC_DIR} -b potential.mac
    DEPENDS potential.mac
)

add_custom_command(OUTPUT pbd_cons.f90
    COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b pbd_constraint.mac
    DEPENDS pbd_constraint.mac
)

add_custom_command(OUTPUT elasticity.f90
    COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b elasticity.mac
    DEPENDS elasticity.mac
)

add_custom_command(OUTPUT stvk_hess_stiff_tensor.f90
    COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b stvk_hess_stiff_tensor.mac
    DEPENDS stvk_hess_stiff_tensor.mac
)

add_custom_command(OUTPUT cosserat_rod.f90
    COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b cosserat_rod.mac
    DEPENDS cosserat_rod.mac
)

add_custom_command(OUTPUT poly_elasticity.f90
    COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b poly_elasticity.mac
    DEPENDS poly_elasticity.mac
)

add_custom_command(OUTPUT nnet.f90
    COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b nnet.mac
    DEPENDS nnet.mac
)

add_custom_command(OUTPUT vox_elastic.f90
    COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b vox_elastic.mac
    DEPENDS vox_elastic.mac
)

add_custom_command(OUTPUT vox_dfdPF.f90
    COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b vox_dfdPF.mac
    DEPENDS vox_dfdPF.mac
)

add_custom_command(OUTPUT vox_elastic_mtr.f90
  COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b vox_elastic_mtr.mac
  DEPENDS vox_elastic_mtr.mac
)

add_custom_command(OUTPUT sdf.f90
  COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b sdf.mac
  DEPENDS sdf.mac
)

add_custom_command(OUTPUT quad_to_lin.f90
  COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b quad_to_lin.mac
  DEPENDS quad_to_lin.mac
)

add_custom_command(OUTPUT fem_basis.f90
  COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b fem_basis.mac
  DEPENDS fem_basis.mac
)

add_custom_command(OUTPUT hex_cons_law.f90
  COMMAND "maxima" ARGS   --userdir=${SRC_DIR} -b hex_cons_law.mac
  DEPENDS hex_cons_law.mac
  )

add_custom_command(OUTPUT coro_basis_def_oper.f90
  COMMAND "maxima" ARGS  --userdir=${SRC_DIR} -b coro_basis_def_oper.mac
  DEPENDS coro_basis_def_oper.mac
  )
# if(CUDA_FOUND)
#     cuda_add_library(cu-jacobi SHARED jacobi_kernel.cu)
#     set(CU_JACOBI_LIBRARY cu-jacobi)
# else(CUDA_FOUND)
#     set(CU_JACOBI_LIBRARY "")
# endif(CUDA_FOUND)

add_library(bigbang SHARED ${source} pbd_cons.f90 elasticity.f90 cosserat_rod.f90 nnls.c poly_elasticity.f90 nnet.f90 vox_elastic.f90 vox_dfdPF.f90 vox_elastic_mtr.f90 sdf.f90 quad_to_lin.f90 fem_basis.f90 hex_cons_law.f90 coro_basis_def_oper.f90)
target_link_libraries(bigbang
    jtf-mesh
    ${SuiteSparse_LIBRARIES}
    ${LBFGS_LIBRARIES}
    ${ARPACK_LIBRARIES}
    ${CU_JACOBI_LIBRARY}
    ann-parallel
    minpack
    hj-math
    ${ALGLIB_LIBRARIES}
    jsoncpp
    zjucad-linear-solver
    ${IPOPT_LIBRARIES}
    jtf-optimizer
    /usr/local/lib/libhsl.so
)
