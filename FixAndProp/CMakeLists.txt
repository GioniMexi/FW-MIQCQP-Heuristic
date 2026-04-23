# Common boilerplate
cmake_minimum_required(VERSION 3.10)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Prevent in-source builds
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory and run CMake from there")
endif()

# Determine if this is built as a subproject (using add_subdirectory)
# or if it is the master project.
set(MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(MASTER_PROJECT ON)
  message(STATUS "Fix and Propagate as master project")
  message(STATUS "CMake version: ${CMAKE_VERSION}")

  project(fix-and-prop VERSION 2.3)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  # # Get the current working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE FIXPROP_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --always --dirty
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE FIXPROP_GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

# Find Threads
find_package(Threads)

# Define libfp (fix and propagate)
add_library(fp SHARED src/transformers.cpp src/c_interface.cpp)
target_link_libraries(fp PUBLIC Utils::Lib fmt::fmt Prop::Lib)
add_library(Fp::Lib ALIAS fp)

target_include_directories(fp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Define fix and prop executable
add_executable(fixprop src/main.cpp)

target_include_directories(fixprop PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

if (APPLE)
	target_link_libraries(fixprop Prop::Lib Fp::Lib Utils::Lib fmt::fmt)
else()
	target_link_libraries(fixprop -Wl,--whole-archive Prop::Lib Fp::Lib -Wl,--no-whole-archive Utils::Lib fmt::fmt)
endif()

# Add subprojects
if (GIT_FOUND AND MASTER_PROJECT)
  add_subdirectory(extern/fmt)
  add_subdirectory(extern/utils)
  add_subdirectory(extern/propagator)
endif()
