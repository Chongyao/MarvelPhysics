file(GLOB header *.h)
file(GLOB source *.cc *.c)

add_library(marval SHARED ${source})

