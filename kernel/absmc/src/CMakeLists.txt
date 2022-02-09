# build script for absmc code.

cmake_minimum_required(VERSION 3.6.3)

project (absmc CXX)

IF (BUILD_MULTISCALE)
    set(CMAKE_PREFIX_PATH "${MUSCLE3_HOME}")
    find_package(PkgConfig REQUIRED)
    pkg_search_module(PKG_MUSCLE3 REQUIRED IMPORTED_TARGET libmuscle)
    #pkg_search_module(PKG_YMMSL REQUIRED IMPORTED_TARGET ymmsl)
ENDIF (BUILD_MULTISCALE)


option (BUILD_APP_3D "Select whether 3D stand-alone applications shall be build." ON)

# set defaults for compiler flags, values are editable by user
set(P_CXX_FLAGS "-fopenmp -std=c++17 -Wall -pedantic -Wno-long-long -Wno-strict-aliasing $ENV{OMP_FLAGS}" CACHE STRING "Flags used by the compiler during all build types.")
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

# make sure our own headers are found
include_directories (.)

#include headers in ISR3D/lib
include_directories (${LOCAL_LIB})

file(GLOB SRC_CORE core/*.h core/*.hh core/*.cpp)
file(GLOB SRC_GRAPHICS graphics/*.h graphics/*.hh graphics/*.cpp)
file(GLOB SRC_KDTREE kdtree/*.h kdtree/*.hh kdtree/*.cpp)
file(GLOB SRC_MODEL_3D model3D/*.h model3D/*.hh model3D/*.cpp)

set (SRC_ABSMC_3D ${SRC_CORE} ${SRC_GRAPHICS} ${SRC_UTIL} ${SRC_KDTREE} ${SRC_MODEL_3D} ${SRC_JAVA} )

#### coast kernel modules ####
IF (BUILD_MULTISCALE)
    message(STATUS "BUILD SMCController")
    add_executable (smc ${SRC_ABSMC_3D} multiscale/SMCController.cpp)
    target_link_libraries (smc PUBLIC PkgConfig::PKG_MUSCLE3)
#    target_link_libraries (smc PUBLIC PkgConfig::PKG_YMMSL)
ENDIF (BUILD_MULTISCALE)

#### stand-alone applications ####

IF (BUILD_APP_3D)
    message(STATUS "BUILD_APP_3D")
    add_library (absmc3D STATIC ${SRC_ABSMC_3D})

    #add_executable (generateCuboid app3D/generateCuboid.cpp)
    #target_link_libraries (generateCuboid absmc3D)

    #add_executable (generateSimpleCylinder app3D/generateSimpleCylinder.cpp)
    #target_link_libraries (generateSimpleCylinder absmc3D)

    add_executable (generateSmoothCylinder app3D/generateSmoothCylinder.cpp)
    target_link_libraries (generateSmoothCylinder absmc3D)

#    add_executable (generateRidgedChannel app3D/generateRidgedChannel.cpp)
#    target_link_libraries (generateRidgedChannel absmc3D)

    add_executable (equilibration3D app3D/equilibration3D.cpp)
    target_link_libraries (equilibration3D absmc3D)

    add_executable (deployment3D app3D/deployment3D.cpp)
    target_link_libraries (deployment3D absmc3D)

    add_executable (deployment3D-stent app3D/deployment3D-stent.cpp)
    target_link_libraries (deployment3D-stent absmc3D)

    add_executable (pressurizationDisplacement3D app3D/pressurizationDisplacement3D.cpp)
    target_link_libraries (pressurizationDisplacement3D absmc3D)

#    add_executable (inflation3D app3D/inflation3D.cpp)
#    target_link_libraries (inflation3D absmc3D)

    add_executable (uniaxialStrainTest app3D/uniaxialStrainTest.cpp)
    target_link_libraries (uniaxialStrainTest absmc3D)

    add_executable (csvArteryLoader3D app3D/csvArteryLoader3D.cpp)
    target_link_libraries (csvArteryLoader3D absmc3D)

 ENDIF (BUILD_APP_3D)
