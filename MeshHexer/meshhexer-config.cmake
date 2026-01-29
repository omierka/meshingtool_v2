
include(CMakeFindDependencyMacro)

find_dependency(OpenMP)
find_dependency(Eigen3 3.4.0)
find_dependency(CGAL 6.0.1)
include(${CMAKE_CURRENT_LIST_DIR}/meshhexerTargets.cmake)
