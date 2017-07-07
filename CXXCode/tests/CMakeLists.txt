configure_file(GPanelMAP.in ${CMAKE_CURRENT_BINARY_DIR}/test_data/GPanelMAP.in COPYONLY)

add_executable(
  GPanelMAP
  GPanelMAP.cpp
  )

target_link_libraries(
  GPanelMAP
  FastPLM
  /usr/local/lib/libarmadillo.dylib
 )
