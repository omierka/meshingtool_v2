# MeshHexer

A geometry grab-bag library.

## Dependencies
To build MeshHexer you need:
- A C++17 compatible C++ compiler
- CMake version 3.28.3 or newer
- [CGAL](www.cgal.org) version 6.0.1 or newer
- [Boost](www.boost.org) version 1.81 or newer
- [PMPLibrary](www.pmp-library.org) version 3.0.0 or newer
 
CGAL, Boost, and PMP will be downloaded and built automatically if they are not found by CMake.

## Building

Assuming you cloned this repo to ``~/meshhexer`` you can build the library via:
```
mkdir build
cd build
cmake ~/meshhexer
cmake --build . --target MeshHexer
```

Additionally run
```
cmake --build . --target meshhexer-cli
```
to build the corresponding command-line application. Or just run ``cmake --build .``.

The following options are supported by the build system (defaults):
- `MESHHEXER_ALLOW_EXTERNAL_DOWNLOAD`: Enable automatic download of missing thirdparty libraries (ON)
- `MESHHEXER_BUILD_APPLICATIONS`: Enable building of applications (ON)
- `MESHHEXER_DOCUMENTATION`: Enable docs target (ON)
- `MESHHEXER_DOCUMENTATION_INTERNAL`: Add internal headers to documentation (OFF)
- `MESHHEXER_HAVE_OMP`: Enable OMP support (ON)
- `MESHHEXER_INSTALL`: Enable installation support (ON)
- `MESHHEXER_PREFER_EXTERNAL_TPL`: Prefer installed thirdparty libraries over downloading missing libraries (ON)
- `MESHHEXER_SHARED_LIB`: Build as a shared library (OFF)
- `MESHHEXER_TESTING`: Enable unit tests (ON)
- `MESHHEXER_TPL_CACHE_DIRECTORY`: Set cache directory for downloaded thirdparty libraries (None)

## Installing

To install the library (and its dependencies, if they were built during the configure step) run
```
cmake --install . --prefix=path/to/install/dir
```

MeshHexer will be installed as a CMake package and can be found via ``find_package(MeshHexer)`` if the install directory is seen by CMake.

## Documentation

MeshHexer is documented using Doxygen. To build the documentation run
```
cmake --build . --target doc
```
in a configured build directory.

By default only documentation for the public API of MeshHexer is build.
You can build the documentation for internal code by settings the `MESHHEXER_DOCUMENTATION_INTERNAL` CMake variable to `ON`.

## Tests

MeshHexer has a unit tests. To run the tests either build the whole project using
```
cmake --build .
```
or just the tests using
```
cmake --build . --target tests
```

Then run `ctest` to run the test suite.

## Using MeshHexer

MeshHexer installs itself as a proper CMake package. After installation it can be found via a ``find_package(MeshHexer)`` call, if it has been installed to the usual system directories or some path in the ``CMAKE_PREFIX_PATH``. 
After the ``find_package`` call you can link against the ``MeshHexer::MeshHexer`` CMake target.

If you have included MeshHexer as part of your build tree, for example by vendoring the library or via CMake's FetchContent mechanism, the ``MeshHexer::MeshHexer`` target will be available as well.

Build the documentation or see ``src/meshhexer.hpp`` for documentation on the public API.

## meshhexer-cli

``meshhexer-cli`` is a command-line application that gives access to MeshHexers's features without having to integrate the library into your own application.

Call as
```
meshhexer-cli [-h|--help] [<global args>] <command> <mesh> [<args>]
```
Running ``meshhexer-cli --help`` gives more detailed information.

The possible commands are
- ``fbm-mesh``: Produce a non-fitted base mesh for a mesh-hierarchy suitable for a multigrid solver.

- ``min-gap``: Heuristically try to find the smallest "real" gap between opposite faces of a surface mesh. Tries to discard gaps caused by local geometry, such as in that sharp edges of the mesh, and gaps caused by badly meshes surfaces.

- ``report``: Generate a report about the given mesh. Includes topological information like whether the mesh is closed and the default orientation of surface normals, as well as defects like degenerate triangles.

- ``warnings``: Generates a list of warnings about bad triangles of the mesh. Currently checks for self-intersecting triangles, degenerate triangles, and highly anisotropic triangles. Pass ``--summarize`` to get a summary of warnings, instead of printing each warning individually.


