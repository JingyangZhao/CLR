cmake_minimum_required(VERSION 3.15)
project(clr)

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME} main.cpp constant.h graph.cpp graph.h tour.cpp tour.h factory.cpp factory.h ts.cpp ts.h ps.cpp ps.h tsp.cpp tsp.h)
