cmake_minimum_required(VERSION 3.6)

project("invdet_aux")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 11)

find_package(loonlib REQUIRED)
if(loonlib_FOUND)
    include_directories(${LOONLIB_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Cannot find 'loonlib' package")
endif()

find_package(bioinfo REQUIRED)
if(bioinfo_FOUND)
    include_directories(${BIOINFO_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Cannot find 'bioinfo' package")
endif()

add_executable(extract_invaln extract_invaln.cpp)
target_link_libraries(extract_invaln ${BIOINFO_LIBRARIES} ${LOONLIB_LIBRARIES})

install(TARGETS extract_invaln DESTINATION bin)


