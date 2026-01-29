find_package(CGAL REQUIRED COMPONENTS Core)

add_library(cgal_wrapper STATIC cgal_interface.cpp)
target_link_libraries(cgal_wrapper PUBLIC CGAL::CGAL CGAL::CGAL_Core)
target_compile_features(cgal_wrapper PRIVATE cxx_std_17)

add_executable(meshcleaner fortran_driver.f90)
target_link_libraries(meshcleaner PRIVATE cgal_wrapper)

install(TARGETS meshcleaner)
