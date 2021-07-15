

add_executable(projnoise main.cpp)


target_include_directories(projnoise PUBLIC ${CMAKE_CURRENT_BINARY_DIR})



target_link_libraries(projnoise PRIVATE polymer)
