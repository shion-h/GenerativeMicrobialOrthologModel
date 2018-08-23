cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)
add_executable(GMOM src/GMOM.cpp src/GibbsSamplerFromGMOM.cpp)
TARGET_LINK_LIBRARIES(GMOM boost_program_options)
