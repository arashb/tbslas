set(example_IDs
  semilag_01
)
set(example_LIBS tricubic)

foreach(example_id ${example_IDs})
  add_executable(${example_id} ${example_id}.cc)
  target_link_libraries(${example_id} ${example_LIBS})
endforeach(example_id)
