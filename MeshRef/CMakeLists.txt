add_executable(meshref
  src/var.f90
  src/inout.f90
  src/def.f90
  src/cleanup.f90
  src/main.f90
  ${CMAKE_SOURCE_DIR}/common/preprocessor_config.f90
)

set_target_properties(meshref
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)

install(TARGETS meshref DESTINATION bin)
install(DIRECTORY VERTEX DESTINATION bin)
