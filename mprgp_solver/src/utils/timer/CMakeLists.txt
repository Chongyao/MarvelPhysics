cmake_minimum_required(VERSION 2.6)
project(timer)

if(CHAOS_UTILS_LIBRARY_FLAG OR CHAOS_TIMER_LIBRARY_FLAG)
  add_definitions(-DCHAOS_TIMER_LIBRARY_FLAG)
  add_library(ctimer SHARED timer_utils.cc)
  add_library(ctimer_static STATIC timer_utils.cc)
endif(CHAOS_UTILS_LIBRARY_FLAG OR CHAOS_TIMER_LIBRARY_FLAG)
