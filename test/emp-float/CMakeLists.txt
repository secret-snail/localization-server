cmake_minimum_required(VERSION 3.12)
set(common_test_dir "${CMAKE_CURRENT_SOURCE_DIR}/../common-test/")
file(GLOB_RECURSE common_test CONFIGURE_DEPENDS "${common_test_dir}/*.cpp" "${common_test_dir}/*.hpp")
list(FILTER common_test EXCLUDE REGEX ".*eth3d_tester\\.cpp$")
list(FILTER common_test EXCLUDE REGEX ".*kitti_tester\\.cpp$")

# required for emp
set(CMAKE_C_FLAGS "-pthread -Wall -march=native -O3 -maes -mrdseed")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

add_executable(emp_float_server_tester
  "${CMAKE_CURRENT_SOURCE_DIR}/test.cpp"
  ${common_test})
target_include_directories(emp_float_server_tester PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
  ${common_test_dir})
target_compile_options(emp_float_server_tester PRIVATE "-Wall" "-Wextra" "-DPPL_GN" "-fsanitize=address")
target_link_options(emp_float_server_tester PRIVATE -fsanitize=address)
target_link_libraries(emp_float_server_tester ${EMP-SH2PC_LIBRARIES} ${OpenCV_LIBS} EmpFloatLocalization Catch2::Catch2WithMain)
# data obliv needs higher timeout
catch_discover_tests(emp_float_server_tester PROPERTIES TIMEOUT 2000)

add_executable(emp_float_eth3d_bench
  "${CMAKE_CURRENT_SOURCE_DIR}/eth3d_bench.cpp"
  ${common_test})
target_include_directories(emp_float_eth3d_bench PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
  ${common_test_dir})
target_compile_options(emp_float_eth3d_bench PRIVATE "-Wall" "-Wextra")
target_link_libraries(emp_float_eth3d_bench ${EMP-SH2PC_LIBRARIES} ${OpenCV_LIBS} EmpFloatLocalization)

add_executable(emp_float_client "${CMAKE_CURRENT_SOURCE_DIR}/client.cpp" ${common_test})
target_include_directories(emp_float_client PUBLIC ${common_test_dir})
target_compile_options(emp_float_client PRIVATE "-Wall" "-Wextra")
target_link_libraries(emp_float_client ${EMP-SH2PC_LIBRARIES} ${OpenCV_LIBS} EmpFloatLocalization Catch2::Catch2WithMain)
