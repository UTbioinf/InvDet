cmake_minimum_required(VERSION 3.6)

project("invdet_proj")

if(SKBUILD)
    message(STATUS "Build Cython module")
    find_package(PythonInterp REQUIRED)
    find_package(PythonLibs REQUIRED)
    find_package(PythonExtensions REQUIRED)
    find_package(Cython REQUIRED)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_subdirectory(source)
