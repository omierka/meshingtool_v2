# Setup thirdparty dependencies
include(FetchContent)
include(FindPackageOverrideHelper)
include(PrintPackageConfig)
include(GetTPL)

set(MAKE_AVAIL_LIST)
mark_as_advanced(MAKE_AVAIL_LIST)

set(FETCHED_PUBLIC_DEPENDENCIES)
mark_as_advanced(FETCHED_PUBLIC_DEPENDENCIES)

get_tpl(
  PACKAGE_NAME Eigen3
  VERSION 3.4.0
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
  URL_HASH MD5=4c527a9171d71a72a9d4186e65bea559
  CONFIG
  PUBLIC
)

get_tpl(
  PACKAGE_NAME Boost
  VERSION 1.88
  URL https://github.com/boostorg/boost/releases/download/boost-1.88.0/boost-1.88.0-cmake.zip
  URL_HASH MD5=419f6a9273cb90d4f3f65ca0ae02cd00
  PUBLIC
)

if(NOT Boost_FOUND)
  # Only build libraries required for CGAL
  #set(BOOST_INCLUDE_LIBRARIES graph heap logic)
  # Make boost create Boost::boost target
  set(BOOST_ENABLE_COMPATIBILITY_TARGETS ON)
endif()

get_tpl(
  PACKAGE_NAME CGAL
  VERSION 6.0.1
  URL https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1.tar.xz
  URL_HASH MD5=944c789615bff14a56d78b398ec2cc49
  SOURCE_SUBDIR Non-Existing
  CONFIG
  PUBLIC
)

if(NOT CGAL_FOUND)
  set(CGAL_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/cgal-src)
  set(CGAL_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/cgal-build)
endif()

if(MESHHEXER_TESTING)
  get_tpl(
    PACKAGE_NAME Catch2
    VERSION 3.8.1
    URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.8.1.zip
    URL_HASH MD5=78dc279759d847f79a6ecd0735e75c58
    CONFIG
    PRIVATE
  )

endif()

if(MAKE_AVAIL_LIST)
  message(STATUS "Fetching the following third-party dependencies:  ${MAKE_AVAIL_LIST}")
  foreach(tpl IN LISTS MAKE_AVAIL_LIST)
    message(STATUS "-------------------------------------------------------------------")
    message(STATUS "- Configuring ${tpl}...")
    message(STATUS "-------------------------------------------------------------------")
    FetchContent_MakeAvailable(${tpl})

    message(STATUS "-------------------------------------------------------------------")
    message(STATUS "- Configuring ${tpl}...done")
    message(STATUS "-------------------------------------------------------------------")
  endforeach()
endif()
