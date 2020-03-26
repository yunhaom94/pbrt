cmake_minimum_required ( VERSION 3.1.0 )
set(ROOT "${CMAKE_CURRENT_LIST_DIR}")
project ( pbrt )

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")

file ( GLOB PBRT_SOURCE
  src/utlis/*
  src/accelerators/*
  src/cameras/*
  src/filters/*
  src/integrators/*
  src/lights/*
  src/materials/*
  src/samplers/*
  src/shapes/*
  src/textures/*
  src/media/*
  )



add_executable( pbrt src/pbrt.cpp  "src/core/ray.cpp" )
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ROOT}/lib/eigen-3.3.7/)
