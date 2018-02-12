cmake_minimum_required(VERSION 3.1.0)

project(FastPLM)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(/usr/local/include)
include_directories(src)

add_library(FastPLM
  src/Common.h
  src/CrushQueue.h
  src/CrushQueue.cpp
  src/FixedEffects.h
  src/FixedEffects.cpp
  src/FixedEffectModel.h
  src/FixedEffectModel.cpp
  src/GPSolver.h
  src/GPSolver.cpp
  src/MatrixReader.h
  src/MatrixReader.cpp
  src/LinearModel.h
  src/LinearModel.cpp
  )

target_compile_definitions(FastPLM PUBLIC BUILD_WITHOUT_R)

add_executable(Debug debug/Debug.cpp)

target_link_libraries(Debug FastPLM /usr/local/lib/libarmadillo.dylib)
