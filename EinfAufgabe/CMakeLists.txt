cmake_minimum_required(VERSION 3.13)
project(EinfAufgabe)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(EinfAufgabe
        graph.cpp
        graph.hpp
        main.cpp
        union_find.cpp
        union_find.hpp tour.cpp tour.hpp tspsolvers.cpp tspsolvers.hpp)
