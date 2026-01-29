cmake_minimum_required(VERSION 3.28.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "MeshHexer build type")
endif()

################
# Setup project
################
option(MESHHEXER_ALLOW_EXTERNAL_DOWNLOAD "enable download of TPLs" ON)
option(MESHHEXER_BUILD_APPLICATIONS "should targets for applications be added?" ON)
option(MESHHEXER_DOCUMENTATION "add doc target for building documentation" ON)
option(MESHHEXER_DOCUMENTATION_INTERNAL "also document internal headers in documentation" OFF)
option(MESHHEXER_HAVE_OMP "Build MeshHexer with OpenMP" ON)
option(MESHHEXER_INSTALL "should installation of MeshHexer be set up?" ON)
option(MESHHEXER_PREFER_EXTERNAL_TPL "try to find external libraries before fetching" ON)
option(MESHHEXER_SHARED_LIB "build as shared library" OFF)
option(MESHHEXER_TESTING "should testing be enabled?" ON)
option(MESHHEXER_TPL_CACHE_DIRECTORY "directory for caching thirdparty dependencies" "")

project(
  MeshHexer
  VERSION 0.1
  DESCRIPTION "Meshing tool for hexahedral meshes"
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Setup thirdparty dependencies.
# Ensures targets Boost::boost and CGAL::CGAL are available to link against
include(Thirdparty)

if(FETCHED_PUBLIC_DEPENDENCIES)
  if(MESHHEXER_INSTALL)
    message(WARNING "Public dependencies ${FETCHED_PUBLIC_DEPENDENCIES} were automatically downloaded. Installation of MeshHexer is disabled because installation of dependencies is not supported!")
  endif()
  set(MESHHEXER_INSTALL OFF CACHE BOOL "should installation of MeshHexer be set up?" FORCE)
endif()

if(MESHHEXER_TESTING)
  enable_testing()
endif()

###################
# Compiler options
###################

# Add library for common compiler options
add_library(meshhexer-settings INTERFACE EXCLUDE_FROM_ALL)

set_target_properties(meshhexer-settings PROPERTIES
  CXX_STANDARD 17
)

target_compile_options(meshhexer-settings INTERFACE
  -Wall
  -Wextra
  -Wundef
  -Wshadow
  -Woverloaded-virtual
  -Wuninitialized
  -Wvla
  -Wdouble-promotion
  -Wformat=2
  -Wnonnull
  -Wswitch-bool
  -Wsizeof-array-argument
  -Wbool-compare
  -Wsuggest-override
  -Wnon-virtual-dtor
  -Wdelete-non-virtual-dtor
  -Wshift-negative-value
  -Wduplicated-cond
  -Wduplicated-branches
  -Wrestrict
  -Wdangling-else
  -Wnonnull
  -Wrestrict
  -Walloc-zero
  -Wparentheses
  -Wpedantic
)

#########
# OpenMP
#########

if(MESHHEXER_HAVE_OMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(meshhexer-settings INTERFACE OpenMP::OpenMP_CXX)
endif()

##########
# Testing
##########

if(MESHHEXER_TESTING)
  # Add Catch2 scripts to module path if integrated via FetchContent
  if(NOT Catch2_FOUND)
    message(STATUS "Adding ${catch2_SOURCE_DIR}/extras to module path")
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
  endif()

  add_custom_target(tests)
endif()

################
# Documentation
################

if(MESHHEXER_DOCUMENTATION)
  find_package(Doxygen)
  set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/doc)
  set(DOXYGEN_EXTRACT_ALL YES)
  set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)

  if(MESHHEXER_DOCUMENTATION_INTERNAL)
    doxygen_add_docs(doc README.md include/meshhexer src)
  else()
    doxygen_add_docs(doc README.md include/meshhexer)
  endif()
endif()


#####################
# Add subdirectories
#####################
add_subdirectory(src)

if(MESHHEXER_BUILD_APPLICATIONS)
  add_subdirectory(applications)
endif()

add_subdirectory("include")

##############
# Config file
##############

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/meshhexer_config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/meshhexer_config.hpp)
target_sources(MeshHexer
  PRIVATE FILE_SET private_headers
  TYPE HEADERS
  BASE_DIRS ${CMAKE_CURRENT_BINARY_DIR}
  FILES ${CMAKE_CURRENT_BINARY_DIR}/meshhexer_config.hpp)

############
# Packaging
############

# Setup alias for inclusion via FetchContent
add_library(MeshHexer::MeshHexer ALIAS MeshHexer)

# Setup install
if(MESHHEXER_INSTALL)
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  # Create package version file
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/meshhexer-config-version.cmake
    COMPATIBILITY AnyNewerVersion
  )

  # Install package version file and package config file
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/meshhexer-config-version.cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/meshhexer-config.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/meshhexer
  )

  # Add target to export
  install(
    TARGETS MeshHexer meshhexer-settings
    EXPORT meshhexerTargets
    # DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/meshhexer
    FILE_SET public_headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  # Generate meshhexerTargets.cmake
  install(
    EXPORT meshhexerTargets
    NAMESPACE MeshHexer::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/meshhexer
  )
endif()
