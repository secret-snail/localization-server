cmake_minimum_required(VERSION 3.12)
file(GLOB emp_float_server_src CONFIGURE_DEPENDS "*.h" "*.cpp")
list(FILTER emp_float_server_src EXCLUDE REGEX ".*server-gn\\.cpp$")
list(FILTER emp_float_server_src EXCLUDE REGEX ".*server-lm\\.cpp$")

# required for emp
set(CMAKE_C_FLAGS "-pthread -Wall -march=native -O3 -maes -mrdseed")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")


add_library(EmpFloatLocalization SHARED ${emp_float_server_src} )
target_include_directories(EmpFloatLocalization PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../common/")
target_compile_options(EmpFloatLocalization PRIVATE "-Wall" "-Wextra" "-fPIC")
target_link_libraries(EmpFloatLocalization ${EMP-SH2PC_LIBRARIES} ${OpenCV_LIBS})


add_executable(lm_emp_float_server "${CMAKE_CURRENT_SOURCE_DIR}/server-lm.cpp" ${emp_float_server_src})
target_include_directories(lm_emp_float_server PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../common/")
target_compile_options(lm_emp_float_server PRIVATE "-Wall" "-Wextra")
target_link_libraries(lm_emp_float_server ${EMP-SH2PC_LIBRARIES} ${OpenCV_LIBS})
