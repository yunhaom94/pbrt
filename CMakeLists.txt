cmake_minimum_required ( VERSION 3.1.0 )

project ( PBRT )

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/")

ADD_EXECUTABLE ( pbrt src/pbrt.cpp )