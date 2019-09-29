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


add_executable(emp_fixed_server_tester "${CMAKE_CURRENT_SOURCE_DIR}/test.cpp" ${common_test})
target_include_directories(emp_fixed_server_tester PUBLIC ${common_test_dir})
target_compile_options(emp_fixed_server_tester PRIVATE "-Wall" "-Wextra")
target_link_libraries(emp_fixed_server_tester ${EMP-SH2PC_LIBRARIES} ${OpenCV_LIBS} EmpFloatLocalization)
