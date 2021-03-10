# build script for collector

cmake_minimum_required(VERSION 3.6.3)

project (collector-cpp CXX)

set(CMAKE_PREFIX_PATH "${MUSCLE3_HOME}")
find_package(PkgConfig REQUIRED)
pkg_search_module(PKG_MUSCLE3 REQUIRED IMPORTED_TARGET libmuscle)

# set defaults for compiler flags, values are editable by user
set(P_CXX_FLAGS "-fopenmp -std=c++14 -Wall -pedantic -Wno-long-long -Wno-strict-aliasing" CACHE STRING "Flags used by the compiler during all build types.")
set(P_CXX_FLAGS_DEBUG "-O -g" CACHE STRING "")
#set(P_CXX_FLAGS_RELEASE "-O3 -funroll-loops -fno-signed-zeros -fno-trapping-math -fassociative-math -DNDEBUG -mtune=native -march=native" CACHE STRING "")
set(P_CXX_FLAGS_RELEASE "-O3 -funroll-loops -DNDEBUG -mtune=native -march=native" CACHE STRING "")
set(P_BUILD_TYPE "Release" CACHE STRING "")

# set actual compiler flags
set(CMAKE_CXX_FLAGS ${P_CXX_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${P_CXX_FLAGS_DEBUG})
set(CMAKE_CXX_FLAGS_RELEASE ${P_CXX_FLAGS_RELEASE})
set(CMAKE_BUILD_TYPE ${P_BUILD_TYPE})

message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message(STATUS "MUSCLE3_HOME=${MUSCLE3_HOME}")
message(STATUS "BUILD_COLLECTOR")

#include headers in ISR3D/lib
include_directories (${LOCAL_LIB})

add_executable (collector ./collector.cpp)
target_link_libraries (collector PUBLIC PkgConfig::PKG_MUSCLE3)

