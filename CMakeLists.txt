cmake_minimum_required (VERSION 2.8.11)
project (myRenderer)
set (srcFiles main.cpp lodepng.cpp)

add_executable(${PROJECT_NAME} ${srcFiles})

## Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    target_link_libraries (${PROJECT_NAME} m)
    set(CMAKE_CXX_FLAGS "-O3 -fopenmp")        ## Optimize
endif()
