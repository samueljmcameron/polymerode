file(GLOB LIB_SOURCES [^.]*.cpp)


add_library(polymer ${LIB_SOURCES})




target_include_directories(polymer INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(polymer Eigen3::Eigen)
