# set(example_IDs
#   semilag_01
#   )
# set(example_LIBS tricubic)

# foreach(example_id ${example_IDs})
#   add_executable(${example_id} ${example_id}.cc)
#   target_link_libraries(${example_id} ${example_LIBS})
# endforeach(example_id)

# include(ExternalProject)
# ExternalProject_Add(semilag_tree
#   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/semilag_tree
#   BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/semilag_tree
#   CONFIGURE_COMMAND ""
#   BUILD_COMMAND make bin/semilag_tree
#   INSTALL_COMMAND mv bin/semilag_tree ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
# )
