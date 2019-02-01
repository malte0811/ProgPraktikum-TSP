cmake_minimum_required(VERSION 3.13)
project(HauptAufgabe)

set(CMAKE_CXX_STANDARD 11)

set(COMMON_SOURCES src/lp/linear_program.cpp include/linear_program.hpp src/lp/branch_and_cut.cpp include/branch_and_cut.hpp include/cut_generator.hpp)

add_executable(BranchBoundTest ${COMMON_SOURCES} src/tests/branch_bound_test.cpp)
target_include_directories(BranchBoundTest PUBLIC include)
target_link_directories(BranchBoundTest PUBLIC libs)
target_link_libraries(BranchBoundTest qsopt lemon)
#TODO somehow have all the settings copy over between targets?
add_executable(BranchCutTest ${COMMON_SOURCES} src/tests/cut_test.cpp)
target_include_directories(BranchCutTest PUBLIC include)
target_link_directories(BranchCutTest PUBLIC libs)
target_link_libraries(BranchCutTest qsopt lemon)