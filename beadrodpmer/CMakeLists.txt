
configure_file(BeadRodPmerConfig.h.in BeadRodPmerConfig.h @ONLY)

add_executable(beadrodpmer-bin main.cpp run.cpp input.cpp)


target_include_directories(beadrodpmer-bin PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(beadrodpmer-bin PRIVATE beadrodpmer)

include(GNUInstallDirs)

install(TARGETS beadrodpmer-bin DESTINATION bin)
