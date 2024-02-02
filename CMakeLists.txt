cmake_minimum_required(VERSION 3.1)
project(MLSD)

find_package(PNG)
if(PNG_FOUND)
  add_definitions(${PNG_DEFINITIONS})
else(PNG_FOUND)
  # zlib is needed by libPNG
  find_package(ZLIB)
  if(NOT ZLIB_FOUND)
    add_subdirectory(third_party/zlib-1.2.12)
    set(ZLIB_LIBRARY ${ZLIB_LIBRARIES})
  endif(NOT ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})

  set(SKIP_INSTALL_ALL 1)
  add_subdirectory(third_party/libpng-1.6.37)
endif(PNG_FOUND)
include_directories(${PNG_INCLUDE_DIRS})

set(LINE_DETECTION
  lsd.c lsd.h
  segment.cpp segment.hpp
  mlsd.cpp mlsd.hpp
  detection.cpp detection.hpp)

# MULTISCALE LSD
add_executable(muLSD ${LINE_DETECTION}
  io_png.c io_png.h
  image.cpp image.h
  cmdLine.h
  main_detection.cpp)

target_compile_options(muLSD PRIVATE -Wall -Wextra)
target_link_libraries(muLSD PRIVATE ${PNG_LIBRARIES})

add_executable(test_pyramid
  test_pyramid.cpp cmdLine.h io_png.h io_png.c image.h image.cpp)
target_compile_options(test_pyramid PRIVATE -Wall -Wextra)
target_link_libraries(test_pyramid PRIVATE ${PNG_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
