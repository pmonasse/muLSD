cmake_minimum_required(VERSION 3.13)
project(MLSD)

find_package(PNG)
if(NOT PNG_FOUND)
  find_package(ZLIB) # zlib is needed by libPNG
  if(NOT ZLIB_FOUND)
        set(ZLIB_BUILD_EXAMPLES OFF)
        add_subdirectory(third_party/zlib-1.3.1)
        set(ZLIB_FOUND TRUE)
        set(PNG_BUILD_ZLIB ON)
        set_target_properties(${ZLIB_LIBRARIES} PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${ZLIB_INCLUDE_DIRS}")
        add_library(ZLIB::ZLIB ALIAS "${ZLIB_LIBRARIES}")
  endif()

  set(PNG_SHARED OFF)
  set(PNG_TESTS OFF)
  set(PNG_TOOLS OFF)
  set(SKIP_INSTALL_ALL 1)
  add_subdirectory(third_party/libpng-1.6.44)
  set_target_properties(${PNG_LIBRARIES} PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${PNG_INCLUDE_DIRS}")
  add_library(PNG::PNG ALIAS "${PNG_LIBRARIES}")
endif()

set(LINE_DETECTION
  lsd.c lsd.h
  segment.cpp segment.hpp
  cluster.cpp cluster.hpp
  mulsd.cpp mulsd.hpp)

# MULTISCALE LSD
add_executable(muLSD ${LINE_DETECTION}
  io_png.c io_png.h
  image.cpp image.h
  cmdLine.h
  main.cpp)

target_link_libraries(muLSD PRIVATE PNG::PNG)

add_executable(test_pyramid
  test_pyramid.cpp cmdLine.h io_png.h io_png.c image.h image.cpp)
target_link_libraries(test_pyramid PRIVATE PNG::PNG)

add_executable(compare_lines compare_lines.cpp)

if(CMAKE_CXX_COMPILER_ID MATCHES "(GNU)|(CLANG)")
  set_target_properties(muLSD PROPERTIES COMPILE_FLAGS "-Wall -Wextra")
  set_target_properties(test_pyramid PROPERTIES COMPILE_FLAGS "-Wall -Wextra")
  set_target_properties(compare_lines PROPERTIES COMPILE_FLAGS "-Wall -Wextra")
endif()
