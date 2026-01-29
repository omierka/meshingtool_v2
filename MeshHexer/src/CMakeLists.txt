set(MESHHEXER_SOURCES
  meshhexer.cpp
  io.cpp
  meshing.cpp
  properties.cpp
)

set(MESHHEXER_PRIVATE_HEADERS
  io.hpp
  meshing.hpp
  properties.hpp
  warnings.hpp
)

set(MESHHEXER_TESTS
  properties-test
  io-test
)

if(MESHHEXER_SHARED_LIB)
  add_library(MeshHexer SHARED)
else()
  add_library(MeshHexer STATIC)
endif()

target_link_libraries(MeshHexer PRIVATE meshhexer-settings)
target_sources(MeshHexer PRIVATE ${MESHHEXER_SOURCES})
target_sources(MeshHexer
  PRIVATE FILE_SET private_headers
  TYPE HEADERS
  FILES ${MESHHEXER_PRIVATE_HEADERS})

find_package(CGAL 6.0.1 CONFIG)
find_package(Eigen3 3.4.0)

target_link_libraries(MeshHexer PRIVATE CGAL::CGAL Eigen3::Eigen)


if(MESHHEXER_TESTING)
  find_package(Catch2)

  include(Catch)
  foreach (test ${MESHHEXER_TESTS})
    add_executable(${test} ${test}.cpp)

    target_link_libraries(${test} PRIVATE MeshHexer Catch2::Catch2WithMain CGAL::CGAL)
    target_include_directories(${test} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    add_dependencies(tests ${test})

    catch_discover_tests(${test})
  endforeach()
endif()
