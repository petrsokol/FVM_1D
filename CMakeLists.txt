cmake_minimum_required(VERSION 3.26)
project(FVM_1D_02_faceFlux)

set(CMAKE_CXX_STANDARD 17)

add_executable(FVM_1D_02_faceFlux main.cpp
        structures/Primitive.h
        structures/Primitive.cpp
        structures/Conservative.h
        structures/Conservative.cpp
        alt.cpp
)
