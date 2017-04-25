project("pyinvdet")

add_cython_target(invdet_core CXX)
add_library(invdet_core MODULE ${invdet_core})
python_extension_module(invdet_core)
target_link_libraries(invdet_core cppcore)

install(TARGETS invdet_core LIBRARY DESTINATION pyinvdet)
