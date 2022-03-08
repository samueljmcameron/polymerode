

add_executable(matrix_definitions main.cpp)


target_include_directories(matrix_definitions PUBLIC ${CMAKE_CURRENT_BINARY_DIR})



target_link_libraries(change_in_constraints PRIVATE polymer)
