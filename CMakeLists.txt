# $Id: CMakeLists.txt 197 2012-02-03 13:04:48Z kato $

cmake_minimum_required(VERSION 2.6)

project(Albastru)

set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

option(USE_PARALLEL "Parallel version" ON)
if(USE_PARALLEL)
  set(CMAKE_C_COMPILER mpicc)
  set(CMAKE_CXX_COMPILER mpicxx)
endif(USE_PARALLEL)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_FLAGS "-Wall")
endif(CMAKE_COMPILER_IS_GNUCC)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall")
endif(CMAKE_COMPILER_IS_GNUCXX)

add_subdirectory(Core)
add_subdirectory(Python)
add_subdirectory(IO)
#add_subdirectory(Comm)
add_subdirectory(Tests)

