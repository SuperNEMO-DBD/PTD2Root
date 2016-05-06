# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 2.8)
project(PTD2Root)

# Modules use Bayeux, so we need to locate this or fail
find_package(Bayeux REQUIRED)

# Ensure our code can see the Bayeux headers
include_directories(${Bayeux_INCLUDE_DIRS})

# Build a dynamic library from our sources
add_library(PTD2Root SHARED ptd2root.h ptd2root.cpp)
