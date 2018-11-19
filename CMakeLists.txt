#cmake_minimum_required(VERSION 3.12)
#project(recomp)

#set(CMAKE_CXX_STANDARD 14)

#include_directories(.)

#add_executable(recomp
#        Makefile
#        recompression.cpp
#        rlslp.cpp
#        rlslp.hpp
#        repair.cpp)

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(recomp)
set(PROJECT_VENDOR "Christopher Osthues")
set(PROJECT_CONTACT "osthues.christopher@web.de")
set(PROJECT_URL "https://gitlab.com/osthues.christopher/recomp.git")
set(PROJECT_DESCRIPTION  "LCE: Longest Common Extensions with Recompression")
set(CMAKE_BUILD_TYPE  "Develop")

add_compile_options(-O3 -W -Wall -Wextra -Wpedantic -Wno-deprecated -ftree-vectorize -march=native)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

include_directories(.)

#add_custom_target(external-downloads)
#cmake_policy(SET CMP0046 OLD)

#include(${CMAKE_SOURCE_DIR}/CMakeSettings.cmake)

# Include header files from directory 'include'
#include_directories(${CMAKE_SOURCE_DIR}/include)

#list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
#include(${CMAKE_SOURCE_DIR}/cmake/Dependency.cmake)

add_executable(recomp2 rlslp.hpp rlslp.cpp recompression.cpp)
target_link_libraries(recomp2 "-lm")
install(TARGETS recomp2 RUNTIME DESTINATION build)
