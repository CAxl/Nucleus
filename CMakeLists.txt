cmake_minimum_required (VERSION 3.21)

project (Nucleus)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory("math")
add_subdirectory("nucleus")

