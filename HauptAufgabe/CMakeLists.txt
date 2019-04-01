cmake_minimum_required(VERSION 3.7)
project(HauptAufgabe)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "-O2")

set(COMMON_SOURCES src/lp/linear_program.cpp include/linear_program.hpp
        src/lp/branch_and_cut.cpp include/branch_and_cut.hpp include/cut_generator.hpp
        src/tsp/subtour_cut_gen.cpp include/subtour_cut_gen.hpp src/tsp/tsp_instance.cpp
        include/tsp_instance.hpp include/lemon_fixes/nagamochi_ibaraki.h
        src/tsp/two_matching_cut_gen.cpp include/two_matching_cut_gen.hpp
        src/tsp/tsp_solvers.cpp include/union_find.hpp src/misc/union_find.cpp
        src/tsp/tsp_solution.cpp include/tsp_solution.hpp include/tsp_solvers.hpp
        src/tsp/tsp_solvers.cpp src/misc/relative_tolerance.cpp include/relative_tolerance.hpp
        src/tsp/comb_cut_gen.cpp include/comb_cut_gen.hpp src/misc/tsp_utils.cpp include/tsp_utils.hpp
        src/tsp/connectivity_cut_gen.cpp include/connectivity_cut_gen.hpp src/tsp/one_tree.cpp include/one_tree.hpp)

include_directories(include)
link_directories(libs)
link_libraries(cplex lemon dl pthread)

add_executable(BranchBoundTest ${COMMON_SOURCES} src/tests/branch_bound_test.cpp)
add_executable(BranchCutTest ${COMMON_SOURCES} src/tests/cut_test.cpp)
add_executable(TSPtest2 ${COMMON_SOURCES} src/tests/tsp_test.cpp)
add_executable(TSPsolver ${COMMON_SOURCES} src/main.cpp)
add_executable(2MatchingSim ${COMMON_SOURCES} src/tests/two_matching_sim.cpp)
add_executable(CombTest ${COMMON_SOURCES} src/tests/comb_test.cpp)
add_executable(1Tree ${COMMON_SOURCES} src/tests/one_tree_test.cpp)
