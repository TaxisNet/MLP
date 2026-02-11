cmake_minimum_required(VERSION 3.10)
project(MLP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(mlp
    src/main.cpp
    src/readData.cpp
    src/search.cpp
    src/construction.cpp
    src/neighborhoods.cpp
)

target_include_directories(mlp PRIVATE src)


# Pybind11 module
option(BUILD_PYBIND "Build pybind11 module" ON)

if(BUILD_PYBIND)
  set(PYBIND11_FINDPYTHON ON)
  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    include(FetchContent)
    FetchContent_Declare(pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG v2.11.1
    )
    FetchContent_MakeAvailable(pybind11)
  endif()

  pybind11_add_module(mlp_py
    src/pybind_module.cpp
    src/search.cpp
    src/construction.cpp
    src/neighborhoods.cpp
  )

  target_include_directories(mlp_py PRIVATE src)
endif()