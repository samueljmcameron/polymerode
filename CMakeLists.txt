cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(beadrodpmer VERSION 1.0 DESCRIPTION "Single tethered polymer dynamics"
        LANGUAGES CXX)


include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)

if( supported )
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
	message(STATUS "IPO / LTO not supported: <${error}>")
endif()


find_package(Eigen3 3.4 REQUIRED NO_MODULE)


add_subdirectory(libbeadrodpmer)
add_subdirectory(beadrodpmer)
#add_subdirectory(tests/noise_projections)
#add_subdirectory(tests/change_in_constraints)
