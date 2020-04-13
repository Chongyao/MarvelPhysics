cmake_minimum_required(VERSION 2.6)
project(tests)

# test mprgp
# add_definitions(-DSINGLE_PRINTER_FLAG)
# add_definitions(-DEMERG_PRINT_FLAG)
# add_definitions(-DALERT_PRINT_FLAG)
# add_definitions(-DCRITICAL_PRINT_FLAG)
# add_definitions(-DERROR_PRINT_FLAG)
# add_definitions(-DWARN_PRINT_FLAG)
# add_definitions(-DNOTICE_PRINT_FLAG)
# add_definitions(-DINFO_PRINT_FLAG)
# add_definitions(-DDEBUG_PRINT_FLAG)
# add_definitions(-DASSERT_EXT_FLAG)

# add_definitions(-DSINGLE_TIMER_FLAG)
# add_definitions(-DUSE_STOPWATCH_HELPER_FLAG)

add_executable(test_mprgp test_mprgp.cc)
target_link_libraries(test_mprgp SQP)

add_executable(test_collid test_collid.cc)
target_link_libraries(test_collid coll_v3)
