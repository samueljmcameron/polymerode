file(GLOB LIB_SOURCES [^.]*.cpp)


add_library(beadrodpmer ${LIB_SOURCES})



target_include_directories(beadrodpmer
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(beadrodpmer Eigen3::Eigen)

set(LIB_HEADERS atom.hpp bond.hpp double_tether.hpp globalparams.hpp input.hpp iovtk.hpp iovtk.tpp
    matrix.hpp no_tether.hpp polymer.hpp single_tether.hpp)

message(STATUS "LIB_HEADERS: ${LIB_HEADERS}")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/includes/beadrodpmer)

foreach(_HEADER ${LIB_HEADERS})
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/includes/beadrodpmer/${_HEADER} COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${_HEADER} ${CMAKE_CURRENT_BINARY_DIR}/includes/beadrodpmer/${_HEADER} DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_HEADER})
    add_custom_target(${_HEADER} DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/includes/beadrodpmer/${_HEADER})
    add_dependencies(beadrodpmer ${_HEADER})
endforeach()
