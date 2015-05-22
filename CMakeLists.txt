cmake_minimum_required (VERSION 2.6)
project (Comp_graphics)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "-g -Wall  -std=c11 -fopenmp" )
endif()



INCLUDE_DIRECTORIES( include )

ADD_EXECUTABLE( render  main.c model.c model.h tga.c tga.h geometry.h geometry.c)

TARGET_LINK_LIBRARIES( render m )

#make cmake work on linux
#sudo apt-get update && sudo apt-get install build-essential