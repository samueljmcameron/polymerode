
configure_file(PolyConfig.h.in PolyConfig.h @ONLY)

add_executable(poly main.cpp)


target_include_directories(poly PUBLIC ${CMAKE_CURRENT_BINARY_DIR})



target_link_libraries(poly PRIVATE polymer)
