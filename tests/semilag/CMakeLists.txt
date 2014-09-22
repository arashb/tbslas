# test variables
# id of each test (.cc file name wihout the .cc extention)
set(test_IDs
  vec_field_test
  traj_test
  utils_test
  cubic_interp_policy_test
  )
set(test_LIBS  gtest gtest_main tricubic)

foreach(test_id ${test_IDs})
  add_executable(${test_id} ${test_id}.cc)
  target_link_libraries(${test_id} ${test_LIBS})
  add_test(${test_id} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_id})
endforeach(test_id)
