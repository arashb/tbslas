##############################
# CMAKE CONFIG
##############################
cmake_minimum_required(VERSION 2.8)

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
include_directories(${tbslad_SOURCE_DIR}/include)
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
# SUB-DIRECTORIES
##############################
add_subdirectory(src)

##############################
# UNITTESTS
##############################
if(BUILD-GTEST-UNITTEST)
  set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
  add_subdirectory(${EXT_PROJECTS_DIR}/gtest)
  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
  add_subdirectory(tests)
endif()