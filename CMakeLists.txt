cmake_minimum_required(VERSION 3.0.0)
project(trace VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "-std=c++17")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_executable(trace trace.cpp trace.H slVector.cpp slVector.H slIO.H Bounds3.hpp BVH.cpp BVH.hpp Ray.H
                Object.H Object.cpp)