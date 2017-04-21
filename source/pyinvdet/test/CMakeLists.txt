cmake_minimum_required(VERSION 3.6.1 FATAL_ERROR)

project("test_rmq")

add_subdirectory(LCA/)

include_directories(${RMQ_lib_SOURCE_DIR})
add_executable(test_rmq test.cpp)
