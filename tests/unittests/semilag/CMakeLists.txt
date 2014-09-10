# test variables
# id of each test (.cpp file name wihout the .cpp extention)
set(test_IDs
  VecFieldTest
  TrajTest
  )
set(test_LIBS  gtest gtest_main tricubic)

foreach(test_id ${test_IDs})
  add_executable(${test_id} ${test_id}.cpp)
  target_link_libraries(${test_id} ${test_LIBS})
  add_test(${test_id} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_id})
endforeach(test_id)
