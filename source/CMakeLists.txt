project("invdet")

add_subdirectory(util/)
add_subdirectory(LCA/)
include_directories(${RMQ_lib_SOURCE_DIR} ${util_lib_SOURCE_DIR})

add_subdirectory(cppcore/)
include_directories(${cppcore_lib_SOURCE_DIR})

add_subdirectory(pysrc/)
