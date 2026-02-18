# CGAL + Fortran Mesh Loader

This project demonstrates how to expose a subset of CGAL’s C++ API to Fortran via `ISO_C_BINDING` and then run a full hex-mesh filtering pipeline in Fortran.  The workflow is split into two parts:

1. `src/cgal_interface.cpp` loads an OFF surface with CGAL, triangulates it, and exposes a small C API to query raw vertex/triangle buffers.
2. `src/fortran_driver.f90` (with helpers in `src/bc_treatment.f90`) reads the hexahedral mesh, calls the CGAL wrapper, filters elements, recomputes boundary flags, classifies boundary faces, and writes all derived outputs.

The driver understands the two supported parametrisations today—**HollowCylinder** and **Box**—and can be extended as new mesh types appear.

## Building

```bash
module load cmake/3.28.3 gcc/latest-v13 openmpi/4.1.6 python/3.13.5 cgal/6.0.1
export MPFR_INCLUDE_DIR=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/include
export MPFR_LIBRARIES=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/lib64/libmpfr.a
export GMP_INCLUDE_DIR=/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/include
export GMP_LIBRARIES="/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/lib64/libgmp.a;/sfw/gcc/13.2.0-static-gmp-mpfr-mpc-isl/lib64/libgmpxx.a"

cmake -S . -B build
cmake --build build
```

The repository ships with `runner.sh`, which loads the same modules (including OpenMPI to provide the Fortran MPI runtime), builds the tree, and runs the default PROFEX example.  Modify the final invocation in `runner.sh` (several examples are commented out) to batch-test different data sets.

## Running

The binary lives at `./build/src/fortran_cgal`.  Useful flags:

| Flag | Meaning |
|------|---------|
| `-h`, `--hex` | Path to the input *.tri hexahedral mesh. |
| `-t`, `--tri` | Path to the OFF surface mesh used for filtering. |
| `-s`, `--hex-scale` | Optional scale applied to the *input* hex mesh prior to filtering (defaults to 1.0). |
| `-o`, `--output-folder` | Directory that contains `MESHCONFIG.dat` and receives outputs.  The driver creates `<output>/meshDir` automatically. |

Example:

```bash
# Parallel run (rank 0 orchestrates, all work happens on ranks >=1)
mpirun -np 8 ./build/src/fortran_cgal \
    -h PROFEX/CRS.tri \
    -t PROFEX/surface.off \
    -s 10.0 \
    -o PROFEX
```

When `-o` is supplied the driver expects `<output>/MESHCONFIG.dat`.  That INI-style file describes the originating mesh (Box or HollowCylinder), the overall bounding box/barrel, and the element spacing.  The parser lives in `bc_treatment.f90` and normalises section/key names, so the same format can be used for future mesh families.

The driver now requires at least **two** MPI ranks: rank 0 handles orchestration/output only, while ranks ≥1 process the hexahedral workload. Rank 0 still reduces the filtered mesh, writes files, and prints the boundary summaries.

## Filtering pipeline

1. The driver reads the hexahedral mesh, rescales the coordinates (if `-s` is set), and loads the OFF surface through CGAL.
2. For each hexahedron a minimal bounding sphere is computed; elements whose spheres intersect the surface, contain vertices inside the surface, or whose centres lie inside the surface are marked as “kept”.
3. The connectivity and coordinates are then reduced to the kept elements. `Filtered.vtu` receives the unscaled coordinates, while `meshDir/Filtered.tri` is written with coordinates multiplied by 0.1 (as required by the consuming tooling).
4. Boundary flags (`KNPR`) are recomputed: every hexahedral face is tracked, and the vertices belonging to faces that only appear once in the mesh receive `KNPR=1`.

## Boundary classification and outputs

The helper module `bc_treatment.f90` owns the classification logic:

- The shortest hexahedral edge in the filtered mesh is located and its half-length becomes the *tolerance*—reported in the console—to decide whether a face belongs to a parametrisation.
- **HollowCylinder**: outer cylinder, inner cylinder, axial min plane, axial max plane.  A face enters a bucket only if *all four* of its vertices satisfy the respective radius/plane equation within the tolerance.  “Inner wall” faces are those boundary faces that don’t match any of the four parametrisations.
- **Box**: x±, y±, z± planes plus inner wall faces.  As with cylinders, all vertices must lie on a plane to classify the face.
- After a face bucket is computed, the corresponding vertex list is derived from those faces—nodal assignments are strictly driven by face membership so the two stay in sync.

### Files emitted into `<output>/meshDir`

* `Filtered.tri` &mdash; filtered hex mesh, coordinates scaled by 0.1.
* `*.par` &mdash; one file per parametrisation.  The first line contains `<count> <keyword>` where the keyword is `Wall`, except for the z+ files which use `Outflow`.  The second line literally contains `" "`, followed by one vertex id per line.
  * Hollow cylinder: `cyl_out.par`, `cyl_in.par`, `z-.par`, `z+.par`, `innerwall.par`.
  * Box: `x-.par`, `x+.par`, `y-.par`, `y+.par`, `z-.par`, `z+.par`, `innerwall.par`.
* `file.prj` &mdash; plain text list referencing the files generated in this run (only the filenames, not full paths).  This is useful for downstream automation.

Outside of `meshDir` we also emit `<output>/Filtered.vtu`, which contains the filtered cells with the original scale and the recomputed `KNPR` array: it’s a convenient ParaView/VTK representation for inspection.

## Summary

Between the CGAL wrapper and the Fortran driver you get a reproducible “load, filter, classify, export” pipeline:

1. Build once via `cmake` (or `runner.sh`).
2. Provide a hexahedral mesh, an OFF surface, and a `MESHCONFIG.dat`.
3. Run `fortran_cgal` with `-h`, `-t`, `-s`, and `-o`.
4. Collect `Filtered.vtu`, `meshDir/Filtered.tri`, all parametrisation `.par` files, and the `file.prj` inventory for downstream tooling.
