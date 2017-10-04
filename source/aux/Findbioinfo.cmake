# - Try to find loonlib
# Once done this will define
# LOONLIB_FOUND - System has loonlib
# LOONLIB_INCLUDE_DIRS - The loonlib include directories
# LOONLIB_LIBRARIES - The libraries need to use loonlib

find_path(BIOINFO_INCLUDE_DIR bioinfo/util.h HINTS ${BIOINFO_ROOT_DIR}/include)

find_library(BIOINFO_LIBRARY NAMES bioinfo HINTS ${BIOINFO_ROOT_DIR}/lib)
find_library(BIOINFO_MUMMER_LIBRARY NAMES bioinfo_mummer HINTS ${BIOINFO_ROOT_DIR}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(bioinfo DEFAULT_MSG BIOINFO_INCLUDE_DIR BIOINFO_LIBRARY BIOINFO_MUMMER_LIBRARY)
mark_as_advanced(BIOINFO_INCLUDE_DIR BIOINFO_LIBRARY BIOINFO_MUMMER_LIBRARY)

set(BIOINFO_INCLUDE_DIRS ${BIOINFO_INCLUDE_DIR})
set(BIOINFO_LIBRARIES ${BIOINFO_LIBRARY} ${BIOINFO_MUMMER_LIBRARY})


