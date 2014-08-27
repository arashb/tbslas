##############################
# CMAKE CONFIG
##############################
cmake_minimum_required(VERSION 2.6)

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

##############################
# SUB-DIRECTORIES
##############################
add_subdirectory(src)
add_subdirectory(tests)