cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(poly VERSION 1.0 DESCRIPTION "Single tethered polymer dynamics"
        LANGUAGES CXX)


add_compile_options(-O3)



find_package(Eigen3 3.3 REQUIRED NO_MODULE)


add_subdirectory(src)
add_subdirectory(main)
add_subdirectory(tests/noise_projections)
add_subdirectory(tests/change_in_constraints)
