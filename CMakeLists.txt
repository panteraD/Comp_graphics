cmake_minimum_required (VERSION 2.6)
project (Comp_graphics)

set(CMAKE_C_FLAGS "-g -Wall  -std=c11")

INCLUDE_DIRECTORIES( include )

ADD_EXECUTABLE( render  main.c model.c model.h tga.c tga.h geometry.h geometry.c)

TARGET_LINK_LIBRARIES( render m )