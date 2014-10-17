##############################
# CMAKE CONFIG
##############################
cmake_minimum_required(VERSION 2.8)
enable_language(CXX)

# check for c++11 flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

##############################
# PROJECT CONFIG
##############################
project(tbslad)
enable_testing()

# compile googletest and defined unittests.
# provides a target to run the unittests with "make test".
option(BUILD-GTEST-UNITTEST "Compile googletest and run the unittests." ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${tbslad_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${tbslad_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${tbslad_SOURCE_DIR}/lib)
set(CMAKE_BUILD_TYPE debug)

# include/lib directories
include_directories(${tbslad_SOURCE_DIR}/include
  ${tbslad_SOURCE_DIR}/src
  $ENV{C_INCLUDE_PATH}
  $ENV{CPLUS_INCLUDE_PATH}
  )
link_directories(${tbslad_SOURCE_DIR}/lib)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

##############################
# EXTERNALS
##############################
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)

# tricubic interpolation lib
add_subdirectory(${EXT_PROJECTS_DIR}/tricubic)
include_directories(${tricubic_SOURCE_DIR})

# google c++ test framework
if(BUILD-GTEST-UNITTEST)
  add_subdirectory(${EXT_PROJECTS_DIR}/gtest)
  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
endif()

##############################
# SUB-DIRECTORIES
##############################
add_subdirectory(src)
add_subdirectory(examples)
if(BUILD-GTEST-UNITTEST)
  add_subdirectory(tests)
endif()
