# Required system nauty/Traces headers and library.
include(FindPackageHandleStandardArgs)
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_search_module(PC_Nauty QUIET nauty libnauty)
endif()
find_path(Nauty_INCLUDE_DIR NAMES traces.h
  HINTS ${PC_Nauty_INCLUDE_DIRS} ENV NAUTYDIR
  PATH_SUFFIXES nauty include include/nauty)
find_path(Nauty_CORE_INCLUDE_DIR NAMES nauty.h
  HINTS "${Nauty_INCLUDE_DIR}" ${PC_Nauty_INCLUDE_DIRS} ENV NAUTYDIR
  PATH_SUFFIXES nauty include include/nauty)
find_library(Nauty_LIBRARY NAMES nauty
  HINTS ${PC_Nauty_LIBRARY_DIRS} ENV NAUTYDIR PATH_SUFFIXES lib lib64)
set(Nauty_VERSION "${PC_Nauty_VERSION}")
if(Nauty_CORE_INCLUDE_DIR)
  file(STRINGS "${Nauty_CORE_INCLUDE_DIR}/nauty.h" _nauty_version_line
    REGEX "^[ \t]*#[ \t]*define[ \t]+NAUTYVERSION[ \t]+\"[0-9]+\\.[0-9]+")
  if(_nauty_version_line)
    string(REGEX MATCH "[0-9]+\\.[0-9]+(\\.[0-9]+)?" Nauty_VERSION "${_nauty_version_line}")
  endif()
endif()
# The public headers may occupy separate directories. Preserve the complete
# advertised search paths and also support installations without pkg-config.
set(Nauty_INCLUDE_DIRS "${Nauty_INCLUDE_DIR}" "${Nauty_CORE_INCLUDE_DIR}" ${PC_Nauty_INCLUDE_DIRS})
list(REMOVE_DUPLICATES Nauty_INCLUDE_DIRS)

# Respect explicit Nauty_INCLUDE_DIR/Nauty_LIBRARY overrides. Use pkg-config
# for auxiliary dependencies, without substituting a different libnauty.
set(_nauty_libraries "${Nauty_LIBRARY}")
foreach(_library IN LISTS PC_Nauty_LINK_LIBRARIES)
  get_filename_component(_name "${_library}" NAME)
  if(NOT _name MATCHES "^(lib)?nauty([.-].*)?$")
    list(APPEND _nauty_libraries "${_library}")
  endif()
endforeach()

if(Nauty_INCLUDE_DIR AND Nauty_CORE_INCLUDE_DIR AND Nauty_LIBRARY)
  cmake_push_check_state(RESET)
  set(CMAKE_REQUIRED_INCLUDES "${Nauty_INCLUDE_DIRS}")
  set(CMAKE_REQUIRED_DEFINITIONS "${PC_Nauty_CFLAGS_OTHER}")
  set(CMAKE_REQUIRED_LIBRARIES "${_nauty_libraries}")
  set(CMAKE_REQUIRED_LINK_OPTIONS "${PC_Nauty_LDFLAGS_OTHER}")
  # Recheck after changing an explicit library or include-directory override.
  file(SHA256 "${Nauty_INCLUDE_DIR}/traces.h" _nauty_header_hash)
  file(SHA256 "${Nauty_CORE_INCLUDE_DIR}/nauty.h" _nauty_core_hash)
  string(SHA256 _nauty_probe_key "${Nauty_INCLUDE_DIRS};${PC_Nauty_CFLAGS_OTHER};${PC_Nauty_LDFLAGS_OTHER};${_nauty_libraries};${Nauty_VERSION};${_nauty_header_hash};${_nauty_core_hash}")
  if(NOT _nauty_probe_key STREQUAL Nauty_PROBE_KEY)
    unset(Nauty_TRACES_LINKS CACHE)
    set(Nauty_PROBE_KEY "${_nauty_probe_key}" CACHE INTERNAL "Nauty probe inputs" FORCE)
  endif()
  check_cxx_source_compiles("
    #include <traces.h>
    int main() {
      sparsegraph graph{};
      DEFAULTOPTIONS_TRACES(options);
      TracesStats stats{};
      Traces(&graph, nullptr, nullptr, nullptr, &options, &stats, nullptr);
      schreier_freedyn();
      return 0;
    }" Nauty_TRACES_LINKS)
  cmake_pop_check_state()
endif()
find_package_handle_standard_args(Nauty
  REQUIRED_VARS Nauty_INCLUDE_DIR Nauty_CORE_INCLUDE_DIR Nauty_LIBRARY Nauty_TRACES_LINKS
  VERSION_VAR Nauty_VERSION
  REASON_FAILURE_MESSAGE
    "Install the nauty/Traces development headers and library, or set CMAKE_PREFIX_PATH to their installation prefix.")
mark_as_advanced(Nauty_INCLUDE_DIR Nauty_CORE_INCLUDE_DIR Nauty_LIBRARY Nauty_TRACES_LINKS)

if(Nauty_FOUND AND NOT TARGET Nauty::Nauty)
  add_library(Nauty::Nauty INTERFACE IMPORTED)
  set_target_properties(Nauty::Nauty PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${Nauty_INCLUDE_DIRS}"
    INTERFACE_COMPILE_OPTIONS "${PC_Nauty_CFLAGS_OTHER}"
    INTERFACE_LINK_LIBRARIES "${_nauty_libraries}"
    INTERFACE_LINK_OPTIONS "${PC_Nauty_LDFLAGS_OTHER}")
endif()
