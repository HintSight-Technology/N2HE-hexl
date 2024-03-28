project(hexl_example LANGUAGES C CXX)
cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 11)

include_directories(include)
include_directories(/usr/include)
include_directories(/usr/local/include)

include(FindPkgConfig)
if(NOT PKG_CONFIG_FOUND)
  message(FATAL_ERROR "pkg-config not found!" )
endif()

pkg_check_modules(HEXL REQUIRED IMPORTED_TARGET hexl)

get_target_property(HEXL_LINK_LIBRARIES PkgConfig::HEXL INTERFACE_LINK_LIBRARIES)
if(${CMAKE_VERSION} VERSION_LESS "3.15.0" AND ${HEXL_LINK_LIBRARIES} MATCHES "debug")
  set_target_properties(PkgConfig::HEXL PROPERTIES INTERFACE_LINK_OPTIONS "-fsanitize=address")
  set_target_properties(PkgConfig::HEXL PROPERTIES INTERFACE_COMPILE_OPTIONS "-fsanitize=address")
endif()

set_target_properties(PkgConfig::HEXL PROPERTIES INTERFACE_COMPILE_OPTIONS "-g")

add_executable(face FacialRecognition.cpp)
target_link_libraries(face pthread crypto PkgConfig::HEXL)

add_executable(mnist MNIST.cpp)
target_link_libraries(mnist pthread crypto PkgConfig::HEXL)

