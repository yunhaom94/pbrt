cmake_minimum_required ( VERSION 3.1.0 )
set(ROOT "${CMAKE_CURRENT_LIST_DIR}")
project ( pbrt )

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/")

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

file ( GLOB PBRT_SOURCE
  include/utlis/*
  include/accelerators/*
  include/cameras/*
  include/filters/*
  include/integrators/*
  include/lights/*
  include/materials/*
  include/samplers/*
  include/shapes/*
  include/textures/*
  include/media/*
  )



add_executable( pbrt src/pbrt.cpp )
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ROOT}/lib/eigen/)
