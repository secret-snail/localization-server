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


find_package(Boost 1.71.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_executable(gn_fixed_server_tester "${CMAKE_CURRENT_SOURCE_DIR}/test.cpp" ${common_test})
target_include_directories(gn_fixed_server_tester PUBLIC ${common_test_dir})
target_compile_options(gn_fixed_server_tester PRIVATE "-Wall" "-Wextra" "-DPPL_GN")
target_link_libraries(gn_fixed_server_tester ${OpenCV_LIBS} EmpFloatLocalization ${Boost_LIBRARIES})

add_executable(lm_fixed_server_tester "${CMAKE_CURRENT_SOURCE_DIR}/test.cpp" ${common_test})
target_include_directories(lm_fixed_server_tester PUBLIC ${common_test_dir})
target_compile_options(lm_fixed_server_tester PRIVATE "-Wall" "-Wextra" "-DPPL_LM")
target_link_libraries(lm_fixed_server_tester ${OpenCV_LIBS} EmpFloatLocalization ${Boost_LIBRARIES})


add_executable(fixed_server_benchmark "${CMAKE_CURRENT_SOURCE_DIR}/benchmark.cpp" ${common_test})
target_include_directories(fixed_server_benchmark PUBLIC ${common_test_dir})
target_compile_options(fixed_server_benchmark PRIVATE "-Wall" "-Wextra" "-DPPL_GN")
target_link_libraries(fixed_server_benchmark ${OpenCV_LIBS} EmpFloatLocalization ${Boost_LIBRARIES})


add_executable(emp_float_vs_fixed_benchmark_add "${CMAKE_CURRENT_SOURCE_DIR}/float_vs_fixed_speed.cpp" ${common_test})
target_include_directories(emp_float_vs_fixed_benchmark_add PUBLIC ${common_test_dir})
target_compile_options(emp_float_vs_fixed_benchmark_add PRIVATE "-Wall" "-Wextra" "-DTEST_ADD")
target_link_libraries(emp_float_vs_fixed_benchmark_add ${OpenCV_LIBS} EmpFloatLocalization ${Boost_LIBRARIES})

add_executable(emp_float_vs_fixed_benchmark_mul "${CMAKE_CURRENT_SOURCE_DIR}/float_vs_fixed_speed.cpp" ${common_test})
target_include_directories(emp_float_vs_fixed_benchmark_mul PUBLIC ${common_test_dir})
target_compile_options(emp_float_vs_fixed_benchmark_mul PRIVATE "-Wall" "-Wextra" "-DTEST_MUL")
target_link_libraries(emp_float_vs_fixed_benchmark_mul ${OpenCV_LIBS} EmpFloatLocalization ${Boost_LIBRARIES})
