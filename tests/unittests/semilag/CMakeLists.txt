add_executable(VecFieldTest VecFieldTest.cpp)
target_link_libraries(VecFieldTest gtest gtest_main tricubic)
add_test(VecFieldTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/VecFieldTest)