cmake_minimum_required(VERSION 3.1)
project(tools)

# Compiler options:
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

include_directories(
  include
  )

add_executable(bctest
  src/bctest.cpp)
