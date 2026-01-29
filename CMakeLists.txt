cmake_minimum_required(VERSION 3.28.0)

project(
  NewGenMeshingTool
  VERSION 1.0
  LANGUAGES C CXX Fortran
)

# Overwrite HexMesher options
set(MESHHEXER_BUILD_APPLICATIONS ON CACHE BOOL "" FORCE)
set(MESHHEXER_DOCUMENTATION OFF CACHE BOOL "" FORCE)
set(MESHHEXER_DOCUMENTATION_INTERNAL OFF CACHE BOOL "" FORCE)
set(MESHHEXER_HAVE_OMP ON CACHE BOOL "" FORCE)
set(MESHHEXER_INSTALL OFF CACHE BOOL "" FORCE)
set(MESHHEXER_PREFER_EXTERNAL_TPL ON CACHE BOOL "" FORCE)
set(MESHHEXER_SHARED_LIB OFF CACHE BOOL "" FORCE) # We want hexmesher-cli to be standalone
set(MESHHEXER_TESTING  OFF CACHE BOOL "" FORCE)

# Add subprojects
# Build MeshHexer first, because MeshHexer will auto-download cgal and boost if required
add_subdirectory(MeshHexer)
add_subdirectory(AxiMesher)
add_subdirectory(MeshCleaner)
add_subdirectory(MonitorGeneratorForMesher)

# Install just meshhexer-cli
install(TARGETS meshhexer-cli DESTINATION bin)

# Install runner
install(PROGRAMS runner.sh DESTINATION bin RENAME newgenmeshingtool)
