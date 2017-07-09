cmake_minimum_required(VERSION 3.1.0)

project(FastPLM)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(/usr/local/include)
include_directories(src)

add_library(FastPLM
  src/Common.h
  src/FixedEffect.h
  src/PlainModel.h
  src/MAP.h
  src/GPSolver.h
  src/GradientDescent.h
  src/MatrixReader.h
  src/FixedEffect.cpp
  src/PlainModel.cpp
  src/MAP.cpp
  src/GPSolver.cpp
  src/GradientDescent.cpp
  src/MatrixReader.cpp
  )

target_compile_definitions(FastPLM PRIVATE BUILD_WITHOUT_R)

# add_executable(
#  Test-BasicGradientDescent
#  tests/BasicGradientDescent.cpp
#  )

# target_link_libraries(
#  Test-BasicGradientDescent
#  FastPLM
#  /usr/local/lib/libarmadillo.dylib
# )

# add_subdirectory(tests)
