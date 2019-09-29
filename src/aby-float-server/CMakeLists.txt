cmake_minimum_required(VERSION 3.12)
file(GLOB aby_float_server_src CONFIGURE_DEPENDS "*.h" "*.cpp")
list(FILTER aby_float_server_src EXCLUDE REGEX ".*server-gn\\.cpp$")
list(FILTER aby_float_server_src EXCLUDE REGEX ".*server-lm\\.cpp$")

add_library(AbyFloatLocalization SHARED ${aby_float_server_src} )
target_include_directories(AbyFloatLocalization PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../common/")
target_compile_options(AbyFloatLocalization PRIVATE "-Wall" "-Wextra" "-fPIC")
target_link_libraries(AbyFloatLocalization ABY::aby ${OpenCV_LIBS})

# TODO - build server bins
