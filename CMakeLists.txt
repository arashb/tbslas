##############################
# CMAKE CONFIG
##############################
cmake_minimum_required(VERSION 2.8)

##############################
# PROJECT CONFIG
##############################
project(tbslad)
enable_testing()

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
add_subdirectory(tests)