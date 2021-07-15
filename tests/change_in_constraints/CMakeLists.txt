

add_executable(change_in_constraints main.cpp)


target_include_directories(change_in_constraints PUBLIC ${CMAKE_CURRENT_BINARY_DIR})



target_link_libraries(change_in_constraints PRIVATE polymer)
