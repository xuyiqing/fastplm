cmake_minimum_required(VERSION 3.1.0)

add_executable(basic-cpp main.cpp)

target_compile_definitions(basic-cpp PUBLIC BUILD_WITHOUT_R)

target_link_libraries(basic-cpp FastPLM /usr/local/lib/libarmadillo.dylib)

copy_csv_file(basic-cpp data.csv)
copy_csv_file(basic-cpp inds.csv)
