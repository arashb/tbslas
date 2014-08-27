add_executable(semilag_tst1 semilag_tst1.cpp)
target_link_libraries(semilag_tst1 semilag)

add_test(semilag_tst1 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/semilag_tst1)