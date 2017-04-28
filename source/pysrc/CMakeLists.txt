project("pyinvdet")

add_cython_target(invdet_core CXX)
add_library(invdet_core MODULE ${invdet_core})
python_extension_module(invdet_core)
target_link_libraries(invdet_core cppcore)

add_cython_target(maxcut CXX)
add_library(maxcut MODULE ${maxcut})
python_extension_module(maxcut)
target_link_libraries(maxcut cppcore)

install(TARGETS invdet_core maxcut LIBRARY DESTINATION source/pyinvdet)
