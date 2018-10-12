
file(GLOB header *.h)
file(GLOB source *.cc *.c)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
message("this is cmkae current source dir ${CMAKE_CURRENT_SOURCE_DIR}  ")
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
