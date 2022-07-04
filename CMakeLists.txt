cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(beadrodpmer VERSION 1.0 DESCRIPTION "Single tethered polymer dynamics"
        LANGUAGES CXX)


add_compile_options(-O3 -Wall -Wextra -g -Wuninitialized)

message("${CMAKE_MODULE_PATH}")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)


add_subdirectory(libbeadrodpmer)
add_subdirectory(beadrodpmer)
#add_subdirectory(tests/noise_projections)
#add_subdirectory(tests/change_in_constraints)
