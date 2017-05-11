cmake_minimum_required(VERSION 3.1.0)

project(FastPLM)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(/usr/local/include)
include_directories(include)

add_library(FastPLM
  include/Common.h
  include/FixedEffect.h
  include/PlainModel.h
  include/MAP.h
  include/GPanelMAP.h
  include/GradientDescent.h
  include/MatrixReader.h
  src/FixedEffect.cpp
  src/PlainModel.cpp
  src/MAP.cpp
  src/GPanelMAP.cpp
  src/GradientDescent.cpp
  src/MatrixReader.cpp
  )

add_executable(
  Test-BasicGradientDescent
  tests/BasicGradientDescent.cpp
  )

target_link_libraries(
  Test-BasicGradientDescent
  FastPLM
  /usr/local/lib/libarmadillo.dylib
 )
