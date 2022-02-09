cmake_minimum_required(VERSION 3.6.3)

project(flow3d)
enable_language(CXX)

set(EXECUTABLE_NAME "FlowController3D")

set(CMAKE_PREFIX_PATH "${MUSCLE3_HOME}")
find_package(PkgConfig REQUIRED)
pkg_search_module(PKG_MUSCLE3 REQUIRED IMPORTED_TARGET libmuscle_mpi)

message(STATUS "CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
endif()
message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")
message(${CMAKE_BUILD_TYPE})

# Compiler flags
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wnon-virtual-dtor")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -DNDEBUG -mtune=native -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")

message("Enabling MPI")
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
    #set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    add_definitions(-DPLB_MPI_PARALLEL)
endif()

option(ENABLE_POSIX "Enable POSIX" ON)

if(ENABLE_POSIX)
    message("Enabling POSIX")
    add_definitions(-DPLB_USE_POSIX)
endif()

###############################################################################
# Palabos Library
###############################################################################

include_directories("${PALABOS_HOME}/src")
include_directories("${PALABOS_HOME}/externalLibraries")
include_directories("${PALABOS_HOME}/externalLibraries/Eigen3")

file(GLOB_RECURSE PALABOS_SRC "${PALABOS_HOME}/src/*.cpp")
file(GLOB_RECURSE EXT_SRC "${PALABOS_HOME}/externalLibraries/tinyxml/*.cpp")

add_library(palabos STATIC ${PALABOS_SRC} ${EXT_SRC})

###############################################################################

#include headers in ISR3D/lib
include_directories (${LOCAL_LIB})

add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp")

# Link with the following libraries
target_link_libraries(${EXECUTABLE_NAME} PUBLIC palabos)
target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${MPI_CXX_LIBRARIES})
target_link_libraries(${EXECUTABLE_NAME} PUBLIC PkgConfig::PKG_MUSCLE3)

