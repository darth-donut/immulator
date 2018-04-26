cmake_minimum_required(VERSION 3.4)
project(immulator)
set(EXE immulator)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(${EXE}
        src/main.cpp
        src/germline.h
        src/germline_factory.cpp
        src/germline_factory.h
        src/germline.cpp src/germline_configuration.cpp
        src/germline_configuration.h src/immutils.h
        src/cxxopts.hpp)

