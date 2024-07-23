# Basic stuff all cmake files should have
cmake_minimum_required(VERSION 3.28)
project(flan VERSION 0.1 LANGUAGES C CXX)

add_compile_options(-Wall -Wextra -Wpedantic -O3)

# Where the source code is all stored
add_subdirectory(src)
