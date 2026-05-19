cmake_minimum_required(VERSION 3.28.0)

project(
  NewGenMeshingTool
  VERSION 1.0
  LANGUAGES C CXX Fortran
)

set(CMAKE_BUILD_TYPE "Release")

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
add_subdirectory(MeshRef)
add_subdirectory(MeshDeform)

# Install just meshhexer-cli
install(TARGETS meshhexer-cli DESTINATION bin)

# Install unified Python driver (replaces runner.sh + run_all_cases.sh)
install(PROGRAMS preprocessor.py DESTINATION bin RENAME gendie_preprocessor)
install(FILES preprocessor.cfg DESTINATION bin)
install(PROGRAMS
  PostProc/render_all_snapshots.py
  PostProc/render_filtered_snapshot.py
  DESTINATION bin
)
install(DIRECTORY EXAMPLE DESTINATION .)

# Install project documentation into the release tree.
install(FILES README.md DESTINATION docs)
install(FILES
  MeshCleaner/README
  MeshHexer/README.md
  DESTINATION docs
)
install(DIRECTORY
  MeshCleaner/docs
  MeshDeform/docs
  MeshRef/docs
  PostProc/docs
  DESTINATION docs
)
install(DIRECTORY MonitorGeneratorForMesher/doc DESTINATION docs)
install(DIRECTORY MeshHexer/docs DESTINATION docs)
