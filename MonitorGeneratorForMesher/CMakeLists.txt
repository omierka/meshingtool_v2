cmake_minimum_required(VERSION 3.18)

project(hex_VS_triangulation_intersection LANGUAGES Fortran CXX)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED COMPONENTS Fortran)

find_package(CGAL 6.0 CONFIG QUIET)

if(NOT CGAL_FOUND)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(CGAL_PKG REQUIRED IMPORTED_TARGET cgal)
  set(CGAL_TARGET PkgConfig::CGAL)
else()
  set(CGAL_TARGET CGAL::CGAL)
endif()

add_library(cgal_interface STATIC cgal_interface.cpp)
target_link_libraries(cgal_interface PUBLIC ${CGAL_TARGET})

add_executable(hex_VS_triangulation_intersection
  main.f90
  inout.f90
  TRIvsTET.f90
  cgal_bindings.f90
)

set_target_properties(hex_VS_triangulation_intersection
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)

target_link_libraries(hex_VS_triangulation_intersection
  PRIVATE
    cgal_interface
    MPI::MPI_Fortran
)

install(TARGETS hex_VS_triangulation_intersection DESTINATION bin)
