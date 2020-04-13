cmake_minimum_required(VERSION 2.6)
project(logger)

if(CHAOS_UTILS_LIBRARY_FLAG OR CHAOS_LOG_LIBRARY_FLAG)
  add_definitions(-DCHAOS_LOG_LIBRARY_FLAG)
  add_library(clog SHARED log_utils.cc)
  add_library(clog_static STATIC log_utils.cc)
endif(CHAOS_UTILS_LIBRARY_FLAG OR CHAOS_LOG_LIBRARY_FLAG)
