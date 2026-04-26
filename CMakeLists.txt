cmake_minimum_required(VERSION 3.16)
project(enhanced_ga CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(generate generate.cpp)
target_compile_options(generate PRIVATE -O2 -Wall)

add_executable(schedule schedule.cpp)
target_compile_options(schedule PRIVATE -O3 -Wall)
