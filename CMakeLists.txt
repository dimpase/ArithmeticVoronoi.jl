cmake_minimum_required(VERSION 3.24)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "Use an out-of-source build: cmake -S . -B build")
endif()
project(ArithmeticVoronoiBackend VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
  #include <format>
  int main() { return std::format(\"{}\", 42) != \"42\"; }
  " AV_HAS_CXX20_FORMAT)
if(NOT AV_HAS_CXX20_FORMAT)
  message(FATAL_ERROR "A C++20 compiler and standard library with std::format are required. On macOS, install current Xcode Command Line Tools or select a newer C++ toolchain.")
endif()

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(GNUInstallDirs)
include(CTest)

find_package(Git REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost CONFIG REQUIRED COMPONENTS serialization)
pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)
find_package(Nauty REQUIRED)

include(PinnedBackend)

# The backend is a template library. Consume its headers directly, without
# configuring its standalone collection of unrelated programs or patching it.
add_library(av_backend_headers INTERFACE)
file(GLOB _av_header_dirs LIST_DIRECTORIES true CONFIGURE_DEPENDS
  "${polyhedral_common_SOURCE_DIR}/src_*"
  "${polyhedral_common_SOURCE_DIR}/basic_common_cpp/src_*")
foreach(_dir IN LISTS _av_header_dirs)
  if(IS_DIRECTORY "${_dir}")
    target_include_directories(av_backend_headers INTERFACE "${_dir}")
  endif()
endforeach()
target_include_directories(av_backend_headers INTERFACE
  "${polyhedral_common_SOURCE_DIR}"
  "${polyhedral_common_SOURCE_DIR}/basic_common_cpp"
  "${polyhedral_common_SOURCE_DIR}/permutalib/src"
  "${polyhedral_common_SOURCE_DIR}/basic_common_cpp/sparse-map/include/tsl"
  "${polyhedral_common_SOURCE_DIR}/basic_common_cpp/robin-map/include/tsl"
  "${polyhedral_common_SOURCE_DIR}/basic_common_cpp/hopscotch-map/include/tsl")
target_compile_features(av_backend_headers INTERFACE cxx_std_20)
target_compile_definitions(av_backend_headers INTERFACE SANITY_CHECK)
target_link_libraries(av_backend_headers INTERFACE
  Eigen3::Eigen Boost::serialization Threads::Threads PkgConfig::GMPXX Nauty::Nauty)

option(AV_ENABLE_SANITIZERS "Instrument the exporter with ASan and UBSan" OFF)
add_executable(AV_ExportPerfectOrbits backend/export_perfect.cpp)
target_link_libraries(AV_ExportPerfectOrbits PRIVATE av_backend_headers)
set_target_properties(AV_ExportPerfectOrbits PROPERTIES
  CXX_EXTENSIONS OFF RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
if(AV_ENABLE_SANITIZERS)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    message(FATAL_ERROR "AV_ENABLE_SANITIZERS requires GCC or Clang")
  endif()
  target_compile_options(AV_ExportPerfectOrbits PRIVATE
    -fsanitize=address,undefined -fno-omit-frame-pointer -UNDEBUG -g1 -Wall -Wextra)
  target_link_options(AV_ExportPerfectOrbits PRIVATE -fsanitize=address,undefined)
endif()

# Optional historical audit executable; never built or installed by default.
add_executable(PERF_SerialPerfectComputation EXCLUDE_FROM_ALL
  "${polyhedral_common_SOURCE_DIR}/src_perfect/PERF_SerialPerfectComputation.cpp")
target_link_libraries(PERF_SerialPerfectComputation PRIVATE av_backend_headers)
set_target_properties(PERF_SerialPerfectComputation PROPERTIES
  CXX_EXTENSIONS OFF RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/backend/export_perfect.cpp")
file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/backend/export_perfect.cpp" AV_EXPORTER_SHA256)
configure_file(cmake/WriteBackendInfo.cmake.in WriteBackendInfo.cmake @ONLY)
configure_file(cmake/InstallBackend.cmake.in InstallBackend.cmake @ONLY)
set_property(TARGET AV_ExportPerfectOrbits APPEND PROPERTY LINK_DEPENDS
  "${CMAKE_CURRENT_BINARY_DIR}/WriteBackendInfo.cmake")
add_custom_command(TARGET AV_ExportPerfectOrbits POST_BUILD
  COMMAND "${CMAKE_COMMAND}"
    "-DAV_EXECUTABLE=$<TARGET_FILE:AV_ExportPerfectOrbits>"
    "-DAV_OUTPUT_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}"
    "-DAV_BUILD_CONFIG=$<CONFIG>"
    -P "${CMAKE_CURRENT_BINARY_DIR}/WriteBackendInfo.cmake"
  VERBATIM)
install(TARGETS AV_ExportPerfectOrbits RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/InstallBackend.cmake")

if(BUILD_TESTING)
  # Python is only a test driver; configuration, compilation and installation
  # require no Python interpreter.
  find_package(Python3 3.8 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    foreach(_case IN ITEMS classical2 D5 D8 D12 D13 D40 rejection)
      add_test(NAME "backend.${_case}"
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test/backend_smoke.py"
          "$<TARGET_FILE:AV_ExportPerfectOrbits>" "${_case}")
      set_tests_properties("backend.${_case}" PROPERTIES RUN_SERIAL TRUE TIMEOUT 180)
    endforeach()
  else()
    message(STATUS "Python 3 not found; exporter smoke tests are unavailable")
  endif()
endif()

message(STATUS "Backend headers: ${polyhedral_common_SOURCE_DIR} (${AV_BACKEND_COMMIT})")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}; config.json and config_fields.json select the same all-field exporter")
