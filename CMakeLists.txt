cmake_minimum_required (VERSION 2.8.11)
project (smallpt)
set (src smallpt.cpp lodepng.cpp)

add_executable(${PROJECT_NAME} ${src})
target_link_libraries (${PROJECT_NAME} m)

