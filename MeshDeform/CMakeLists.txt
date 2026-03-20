find_package(CGAL REQUIRED COMPONENTS Core)
find_package(MPI REQUIRED)

add_executable(meshdeform
  src/def.f90
  src/inout.f90
  src/main.f90
  ${CMAKE_SOURCE_DIR}/common/setupe3dfile_reader.f90
  ${CMAKE_SOURCE_DIR}/common/preprocessor_config.f90
)

set_target_properties(meshdeform
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)

target_link_libraries(meshdeform PRIVATE CGAL::CGAL CGAL::CGAL_Core cgal_wrapper)
target_link_libraries(meshdeform PRIVATE MPI::MPI_Fortran)

install(TARGETS meshdeform DESTINATION bin)
