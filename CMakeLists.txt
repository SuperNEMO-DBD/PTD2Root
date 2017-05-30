# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(PTD2Root VERSION 1.0.0)

# Will install module and support program
include(GNUInstallDirs)

# Module uses Falaise, so we need to locate this or fail
find_package(Falaise REQUIRED)

# Build a dynamic library from our sources
add_library(PTD2Root SHARED ptd2root.h ptd2root.cpp)

# Link it to the FalaiseModule target to pick up dependencies
target_link_libraries(PTD2Root Falaise::FalaiseModule)

# - Don't link directly as Bayeux/Falaise will export the
#   needed symbols.
# Ensure dynamic lookup on OS X
if(APPLE)
  set_target_properties(PTD2Root
    PROPERTIES LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()

# - configure the example pipeline file
configure_file(p2r.conf p2r.conf @ONLY)

install(TARGETS PTD2Root DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(PROGRAMS ${PROJECT_SOURCE_DIR}/flptd2root.py DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES ${PROJECT_BINARY_DIR}/p2r.conf DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES readme.pdf DESTINATION ${CMAKE_INSTALL_DOCDIR})
