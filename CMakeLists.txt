cmake_minimum_required(VERSION 3.1.0)

project(FastPLM)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(/usr/local/include)
include_directories(src)

add_library(FastPLM
  src/Common.h
  src/CrushQueue.h
  src/FixedEffects.h
  src/FixedEffectModel.h
  src/GPSolver.h
  src/MatrixReader.h
  src/LinearModel.h
  src/CrushQueue.cpp
  src/FixedEffects.cpp
  src/FixedEffectModel.cpp
  src/GPSolver.cpp
  src/MatrixReader.cpp
  src/LinearModel.cpp
  )

target_compile_definitions(FastPLM PUBLIC BUILD_WITHOUT_R)

add_executable(Debug debug/Debug.cpp)

target_link_libraries(Debug FastPLM /usr/local/lib/libarmadillo.dylib)
