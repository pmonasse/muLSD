cmake_minimum_required(VERSION 3.1)
PROJECT(MLSD)

#FIND_PACKAGE(OpenCV REQUIRED)
OPTION (USE_OpenMP "Use OpenMP" ON)
if(USE_OpenMP)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
endif()
INCLUDE_DIRECTORIES(. DETECTION)
		
SET(LINE_DETECTION
  DETECTION/lsd.cpp	DETECTION/lsd.hpp
  DETECTION/segment.hpp DETECTION/segment.cpp
  DETECTION/mlsd.cpp	DETECTION/mlsd.hpp
  DETECTION/detection.cpp DETECTION/detection.hpp
)
	
# MULTISCALE LSD
ADD_EXECUTABLE(LSD_MULTISCALE ${LINE_DETECTION}
		           main_detection.cpp
		           lodepng/lodepng.cpp)
#TARGET_LINK_LIBRARIES( ${OpenCV_LIBS})

target_compile_options(LSD_MULTISCALE PRIVATE -Wall ${OpenMP_CXX_FLAGS})
target_link_libraries(LSD_MULTISCALE PRIVATE ${OpenMP_CXX_FLAGS})

SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
