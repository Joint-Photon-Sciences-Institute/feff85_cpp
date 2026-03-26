cmake_minimum_required(VERSION 3.20)
project(feff85exafs VERSION 8.5.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler settings
include(cmake/CompilerSettings.cmake)

# Third-party dependencies
add_subdirectory(third_party)

# Core include directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Libraries (in dependency order)
add_subdirectory(src/math)
add_subdirectory(src/par)
add_subdirectory(src/common)
add_subdirectory(src/atom)
add_subdirectory(src/fovrg)
add_subdirectory(src/exch)
add_subdirectory(src/json_io)
add_subdirectory(src/rdinp)
add_subdirectory(src/pot)
add_subdirectory(src/fms)
add_subdirectory(src/xsph)
add_subdirectory(src/path)
add_subdirectory(src/genfmt)
add_subdirectory(src/ff2x)
add_subdirectory(src/debye)
add_subdirectory(src/pipeline)

# Standalone executables
add_subdirectory(apps)

# Tests
enable_testing()
add_subdirectory(tests)
