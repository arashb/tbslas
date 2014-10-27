project(tricubic CXX C)
cmake_minimum_required(VERSION 2.6)

set(tricubic_SRCS
  TriCubicInterpolator.cc
)

# find_package(Boost)
# if(Boost_FOUND)
# include_directories(${Boost_INCLUDE_DIRS})
add_library(tricubic ${tricubic_SRCS})
# target_link_libraries(tricubic ${Boost_LIBRARIES})
# endif()