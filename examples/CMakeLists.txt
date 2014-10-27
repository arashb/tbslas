# set(example_IDs
#   semilag_01
#   )
# set(example_LIBS tricubic)

# foreach(example_id ${example_IDs})
#   add_executable(${example_id} ${example_id}.cc)
#   target_link_libraries(${example_id} ${example_LIBS})
# endforeach(example_id)

# find_package(MPI)
# if (MPI_FOUND)
#   include_directories(${MPI_INCLUDE_PATH})
#   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_FLAGS}")
#   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
# endif()


# find_package(OpenMP)
# if (OPENMP_FOUND)
#   include_directories(${OpenMP_INCLUDE_PATH})
#   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# endif()

# set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# include(ExternalProject)
# ExternalProject_Add(semilag_tree
#   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/semilag_tree
#   BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/semilag_tree
#   CONFIGURE_COMMAND ""
#   BUILD_COMMAND make bin/semilag_tree
#   INSTALL_COMMAND mv bin/semilag_tree ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
# )
