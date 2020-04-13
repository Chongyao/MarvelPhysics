file(GLOB header *.h)
file(GLOB source *.cpp)

add_library(SQP STATIC ${source})
