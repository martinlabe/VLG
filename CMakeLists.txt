cmake_minimum_required(VERSION 3.19)
project(VLG C)

set(CMAKE_C_STANDARD 99)
find_package(igraph REQUIRED)
add_executable(VLG main.c)
target_link_libraries(VLG PUBLIC igraph::igraph)