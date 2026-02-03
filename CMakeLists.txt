cmake_minimum_required(VERSION 3.10)
project(MLP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(mlp
    src/main.cpp
    src/readData.cpp
)

target_include_directories(mlp PRIVATE src)