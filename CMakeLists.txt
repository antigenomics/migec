cmake_minimum_required(VERSION 3.20)
project(migec LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
if(NOT MSVC)
  # MSVC Release is already /O2 and rejects -O3; we do not ship Windows wheels, but a source
  # build on Windows should still configure.
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

option(MIGEC_TESTS  "Build C++ tests"           OFF)
option(MIGEC_BENCH  "Build C++ benchmarks"      OFF)
option(MIGEC_PYTHON "Build the pybind11 module" OFF)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)  # gzipped FASTQ, and the .mig block codec

add_library(migec_core STATIC
  src/types.cpp
  src/fastq.cpp
  src/mig_record.cpp
)
target_include_directories(migec_core PUBLIC include PRIVATE src)
target_link_libraries(migec_core PUBLIC Threads::Threads ZLIB::ZLIB)

if(MIGEC_PYTHON)
  set(PYBIND11_FINDPYTHON ON)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_core src/_bindings.cpp)
  target_link_libraries(_core PRIVATE migec_core)
  install(TARGETS _core DESTINATION migec)
endif()

if(MIGEC_TESTS)
  enable_testing()
  add_executable(migec_tests
    tests/cpp/test_types.cpp
    tests/cpp/test_fastq.cpp
    tests/cpp/test_mig_record.cpp
  )
  target_link_libraries(migec_tests PRIVATE migec_core)
  target_include_directories(migec_tests PRIVATE tests/cpp)
  add_test(NAME migec_tests COMMAND migec_tests)
endif()
