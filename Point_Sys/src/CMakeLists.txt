file(GLOB header *.h)
file(GLOB source *.cc *.c)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
message("this is cmkae current source dir ${CMAKE_CURRENT_SOURCE_DIR}  ")

add_library(marval SHARED ${source})
target_link_libraries(marval
    zjucad-matrix
)

