if(NOT TARGET CGAL::CGAL)
  FetchContent_GetProperties(cgal)
  find_package(Boost 1.81 REQUIRED)
  find_package(Eigen3 3.4.0 CONFIG REQUIRED)

  add_library(meshhexer-cgal-extern INTERFACE)
  target_include_directories(meshhexer-cgal-extern INTERFACE "${cgal_SOURCE_DIR}/include")
  target_link_libraries(meshhexer-cgal-extern INTERFACE Boost::boost Eigen3::Eigen)
  target_compile_definitions(meshhexer-cgal-extern INTERFACE CGAL_DISABLE_GMP CGAL_EIGEN3_ENABLED)
  target_compile_options(meshhexer-cgal-extern INTERFACE -w)

  add_library(CGAL::CGAL ALIAS meshhexer-cgal-extern)
endif()
