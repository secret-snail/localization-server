cmake_minimum_required(VERSION 3.12)
set(common_test_dir "${CMAKE_CURRENT_SOURCE_DIR}/")
file(GLOB_RECURSE common_test CONFIGURE_DEPENDS "${common_test_dir}/*.cpp" "${common_test_dir}/*.h")
list(FILTER common_test EXCLUDE REGEX ".*eth3d_tester\\.cpp$")
list(FILTER common_test EXCLUDE REGEX ".*kitti_tester\\.cpp$")

add_executable(eth3d_tester "${CMAKE_CURRENT_SOURCE_DIR}/eth3d_tester.cpp" ${common_test})
target_include_directories(eth3d_tester PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(eth3d_tester PRIVATE "-Wall" "-Wextra")
target_link_libraries(eth3d_tester ${OpenCV_LIBS})


# WARNING - NOT IMPLEMENTED
#add_executable(kitti_tester "${CMAKE_CURRENT_SOURCE_DIR}/kitti_tester.cpp")
#target_include_directories(kitti_tester PUBLIC)
#target_compile_options(kitti_tester PRIVATE "-Wall" "-Wextra")
#target_link_libraries(kitti_tester ${OpenCV_LIBS})

