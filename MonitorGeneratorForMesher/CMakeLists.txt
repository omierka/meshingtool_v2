find_package(MPI REQUIRED)

add_executable(hex_VS_triangulation_intersection main.f90 inout.f90 TRIvsTET.f90)
set_target_properties(hex_VS_triangulation_intersection PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
target_link_libraries(hex_VS_triangulation_intersection PRIVATE MPI::MPI_Fortran)


install(TARGETS hex_VS_triangulation_intersection DESTINATION bin)
