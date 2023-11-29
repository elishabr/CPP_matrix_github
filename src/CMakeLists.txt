cmake_minimum_required(VERSION 3.0)

project("s21_Matrix_CPP")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MY_SOURCE1 "test_folder/test1.cpp")
set(MY_SOURCE2 "My_s21_matrix_oop/s21_matrix_oop.cpp")


enable_testing()
add_subdirectory(googletest)
add_executable(s21_test 
    ${MY_SOURCE1}
    ${MY_SOURCE2}
)
target_link_libraries(s21_test gtest_main gmock_main)