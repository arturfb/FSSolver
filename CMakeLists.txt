cmake_minimum_required(VERSION 2.8)
include (InstallRequiredSystemLibraries)
#set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
include(CPack)

project(aac-fs)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(CMakeBuildTypes)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-pragmas -fPIC -DIL_STD -fexceptions -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")

find_package(Boost 1.50.0 COMPONENTS program_options timer system chrono filesystem)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
#include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})

add_executable(aac-fs Search.cpp)
target_link_libraries(aac-fs ${Boost_LIBRARIES} rt)