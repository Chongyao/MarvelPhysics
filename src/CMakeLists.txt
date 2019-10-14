file(GLOB header *.h)
file(GLOB source *.cc *.c *.cpp *.cu)

# set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
# message("this is cmkae current source dir ${CMAKE_CURRENT_SOURCE_DIR}  ")
# if(CUDA_FOUND)
#     cuda_add_library(cu-jacobi SHARED jacobi_kernel.cu)
#     set(CU_JACOBI_LIBRARY cu-jacobi)
# else(CUDA_FOUND)
#     set(CU_JACOBI_LIBRARY "")
# endif(CUDA_FOUND)

if (${APPLE})
    add_library(marvel-core SHARED ${source})
    target_link_libraries(marvel-core	
        # coll_rigid
    )
else ()
    cuda_add_library(marvel-core SHARED ${source})
    target_link_libraries(marvel-core	
      coll_rigid
      ${OpenMP_TARGET}


    )
endif(${APPLE})

