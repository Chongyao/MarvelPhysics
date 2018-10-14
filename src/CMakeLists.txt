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

add_library(common SHARED ${source})
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
