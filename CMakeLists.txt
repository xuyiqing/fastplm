cmake_minimum_required(VERSION 3.1.0)

project(FastPLM)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(/usr/local/include)
include_directories(src)

add_library(FastPLM
  src/Common.h
  src/CrashQueue.h
  src/FixedEffect.h
  src/PlainModel.h
  src/FastFESolver.h
  src/GPSolver.h
  src/SlowFESolver.h
  src/MatrixReader.h
  src/FixedEffect.cpp
  src/PlainModel.cpp
  src/FastFESolver.cpp
  src/GPSolver.cpp
  src/SlowFESolver.cpp
  src/MatrixReader.cpp
  )

target_compile_definitions(FastPLM PUBLIC BUILD_WITHOUT_R)

add_executable(Debug debug/Debug.cpp)

target_link_libraries(Debug FastPLM /usr/local/lib/libarmadillo.dylib)
