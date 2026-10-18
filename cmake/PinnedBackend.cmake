include(FetchContent)

set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/backend/pins.json")
file(READ "${PROJECT_SOURCE_DIR}/backend/pins.json" AV_BACKEND_PINS)
string(JSON AV_BACKEND_REPOSITORY GET "${AV_BACKEND_PINS}" repository)
string(JSON AV_BACKEND_COMMIT GET "${AV_BACKEND_PINS}" commit)

FetchContent_Declare(polyhedral_common
  GIT_REPOSITORY "${AV_BACKEND_REPOSITORY}"
  GIT_TAG "${AV_BACKEND_COMMIT}"
  GIT_SUBMODULES_RECURSE TRUE
  UPDATE_DISCONNECTED TRUE
  # Only headers are used. This directory intentionally has no CMake project;
  # FetchContent therefore does not add the upstream standalone build.
  SOURCE_SUBDIR _arithmetic_voronoi_headers_only)
FetchContent_MakeAvailable(polyhedral_common)

function(av_check_revision directory expected)
  execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
    WORKING_DIRECTORY "${directory}"
    RESULT_VARIABLE _status OUTPUT_VARIABLE _actual OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT _status EQUAL 0 OR NOT _actual STREQUAL expected)
    message(FATAL_ERROR "Unpinned backend source at ${directory}: expected ${expected}, got ${_actual}")
  endif()
  execute_process(COMMAND "${GIT_EXECUTABLE}" diff --quiet HEAD --
    WORKING_DIRECTORY "${directory}" RESULT_VARIABLE _dirty)
  if(NOT _dirty EQUAL 0)
    message(FATAL_ERROR "Backend source has tracked modifications: ${directory}. Use a clean pinned checkout.")
  endif()
endfunction()

av_check_revision("${polyhedral_common_SOURCE_DIR}" "${AV_BACKEND_COMMIT}")
string(JSON _count LENGTH "${AV_BACKEND_PINS}" submodules)
math(EXPR _last "${_count} - 1")
foreach(_index RANGE 0 ${_last})
  string(JSON _path MEMBER "${AV_BACKEND_PINS}" submodules ${_index})
  string(JSON _commit GET "${AV_BACKEND_PINS}" submodules "${_path}")
  av_check_revision("${polyhedral_common_SOURCE_DIR}/${_path}" "${_commit}")
endforeach()
