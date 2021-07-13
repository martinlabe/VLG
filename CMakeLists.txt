cmake_minimum_required(VERSION 3.19)
project(VLG C)

set(CMAKE_C_STANDARD 99)
find_package(igraph REQUIRED)

set(SRCS
    main.c
    algo.c
    IO.c
    utils.c
    )

set(HEADERS
    algo.h
    IO.h
    utils.h
    )



add_executable(VLG  ${SRCS}  ${HEADERS})
target_link_libraries(VLG PUBLIC igraph::igraph)
