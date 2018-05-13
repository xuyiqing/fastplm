cmake_minimum_required(VERSION 3.1.0)

function(copy_csv_file target csv)
  add_custom_command(
  TARGET ${target} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
          ${CMAKE_SOURCE_DIR}/tests/manual/${target}/${csv}
          $<TARGET_FILE_DIR:${target}>/${csv})
endfunction(copy_csv_file)

add_subdirectory(basic-cpp)
