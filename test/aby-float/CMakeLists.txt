cmake_minimum_required(VERSION 3.12)
set(common_test_dir "${CMAKE_CURRENT_SOURCE_DIR}/../common-test/")
file(GLOB_RECURSE common_test CONFIGURE_DEPENDS "${common_test_dir}/*.cpp" "${common_test_dir}/*.hpp")
list(FILTER common_test EXCLUDE REGEX ".*eth3d_tester\\.cpp$")
list(FILTER common_test EXCLUDE REGEX ".*kitti_tester\\.cpp$")


add_executable(aby_float_test
  "${CMAKE_CURRENT_SOURCE_DIR}/test.cpp"
  ${common_test})
target_include_directories(aby_float_test PUBLIC
  ${common_test_dir}
  "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(aby_float_test PRIVATE "-Wall" "-Wextra")
target_link_libraries(aby_float_test ABY::aby ${OpenCV_LIBS} AbyFloatLocalization Catch2::Catch2WithMain)
catch_discover_tests(aby_float_test PROPERTIES TIMEOUT 4000)


add_executable(aby_float_eth3d_bench
  "${CMAKE_CURRENT_SOURCE_DIR}/eth3d_bench.cpp"
  ${common_test})
target_include_directories(aby_float_eth3d_bench PUBLIC
  ${common_test_dir}
  "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(aby_float_eth3d_bench PRIVATE "-Wall" "-Wextra")
target_link_libraries(aby_float_eth3d_bench ABY::aby ${OpenCV_LIBS} AbyFloatLocalization)