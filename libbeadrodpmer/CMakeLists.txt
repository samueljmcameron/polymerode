file(GLOB LIB_SOURCES [^.]*.cpp)


add_library(beadrodpmer STATIC ${LIB_SOURCES})



target_include_directories(beadrodpmer
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(beadrodpmer PUBLIC Eigen3::Eigen)


set_target_properties(beadrodpmer PROPERTIES PUBLIC_HEADER "double_tether.hpp;no_tether.hpp;polymer.hpp;single_tether.hpp;initialise.hpp;iovtk.hpp")


include(GNUInstallDirs)

install(TARGETS beadrodpmer LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/beadrodpmer)


# generate pc file for pkg-config

set(target1 beadrodpmer)
configure_file(beadrodpmer.pc.in ${CMAKE_CURRENT_BINARY_DIR}/pkgconfig/beadrodpmer.pc @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pkgconfig/beadrodpmer.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)
