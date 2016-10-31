cmake_minimum_required(VERSION 3.6)
project(CFD_HW)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp BoundaryCond.cpp BoundaryCond.h IntCond.cpp IntCond.h RuntimeParameters.cpp RuntimeParameters.h RuntimeJoke.cpp RuntimeJoke.h SolutionProcedure.cpp SolutionProcedure.h)
add_executable(CFD_HW ${SOURCE_FILES})