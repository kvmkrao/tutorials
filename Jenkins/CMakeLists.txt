cmake_minimum_required(VERSION 3.0)
project(jenkinsexample)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)

add_executable(jenkinsexample main.cpp helloWorld.cpp helloWorld.h)

include(CTest)
add_subdirectory(test)

