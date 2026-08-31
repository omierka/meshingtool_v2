# New Gen Meshing Tool

## Dependencies
Required:
- C++ and Fortran compilers
- CMake (3.28.0 or greater)
- MPI compiler wrappers and runtime (`mpicc`, `mpicxx`, `mpifort`, `mpirun`)

Optional, will be downloaded if missing:
- CGAL (6.0.0 or greater)
- Boost (1.88 or greater)

## Building

Extract the NewGenMeshingTool archive.
Create and enter a build directory:

```
mkdir path/to/build/directory && cd path/to/build/directory
```

Call CMake to set up the build system:

```
cmake path/to/root/of/archive
```

Then build the project via:
```
cmake --build .
```

### Verified Local Build

On the LSIII cluster environment, use the module stack below. It provides a
new enough CMake, matching GCC/OpenMPI compiler wrappers, and the CGAL/Boost
packages used by the current build.

Run these commands from the repository root:

```
source /etc/profile.d/modules.sh
module purge
module load cmake/3.28.3 gcc/latest-v13 openmpi/4.1.6 python/3.13.5 cgal/6.0.1 boost/1.88

cmake -S . -B build-ngmt \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpifort

cmake --build build-ngmt -j 8
```

The explicit MPI compiler wrappers are recommended because several components
are MPI programs or link against MPI-enabled Fortran code.

## Installation

Run
```
cmake --install .
```
to install the tool to the default location.

Run
```
cmake --install . --prefix=/path/to/install/location
```
to install the tool to a custom location. Then add the installed `bin/`
directory to your `PATH` via
```
export PATH=/path/to/install/location/bin:$PATH
```

For a local install below the repository root, run:

```
cmake --install build-ngmt --prefix install-ngmt
```

The installed `bin/` directory contains the driver and the runtime tools it
launches:

- `gendie_preprocessor`
- `meshhexer-cli`
- `generate_hollow_cylinder_mesh`
- `meshcleaner`
- `hex_VS_triangulation_intersection`
- `meshref`
- `meshdeform`

The Python driver resolves these tools relative to its own installed `bin/`
directory, so the most reliable way to run installed examples is from the
install root.

## Usage
The main end-to-end driver is `gendie_preprocessor` after installation, or
`preprocessor.py` when running from the source tree.

For a practical end-to-end example from the source tree, run:
```
python3 preprocessor.py case -f EXAMPLE -n <num_proc>
```

This uses the bundled [`EXAMPLE`](EXAMPLE) case and is the recommended first run for new users.
For cluster-oriented runtime guidance, see [`docs/SingleNodeClusterGuide.md`](docs/SingleNodeClusterGuide.md).

An input case folder for the preprocessing workflow needs to contain at least:
- ```surface.off```, the surface mesh to generate a coarse mesh for
- ```setup.e3d```, preprocessing and simulation settings

### First Practical Run

The repository ships with a ready-to-use example folder:

- [`EXAMPLE/setup.e3d`](EXAMPLE/setup.e3d)
- [`EXAMPLE/surface.off`](EXAMPLE/surface.off)

After building and installing into `install-ngmt`, run the installed example
from the install root:

```
cd install-ngmt

source /etc/profile.d/modules.sh
module purge
module load cmake/3.28.3 gcc/latest-v13 openmpi/4.1.6 python/3.13.5 cgal/6.0.1 boost/1.88

export OMP_NUM_THREADS=2
export OMP_PROC_BIND=close
export OMP_PLACES=cores

bin/gendie_preprocessor case -f EXAMPLE -n 2 --skip-modules
```

This runs the full preprocessing workflow on the bundled example case. The
`-n 2` option uses two MPI ranks; this is the minimum practical value because
`meshcleaner` is MPI-based. The `--skip-modules` option tells the driver not to
reload modules internally because the environment was already prepared above.

By default, MPI stages use `mpirun -np <num_proc>`. Within an existing Slurm
allocation, pass `--use-srun` to launch them as `srun <executable>` instead:

```
bin/gendie_preprocessor case --use-srun -f EXAMPLE -n 2 --skip-modules
```

In `--use-srun` mode, the task count is inherited from the Slurm allocation;
`-n` is not added to the `srun` command.

For a larger single-node run, increase both the OpenMP thread count and the MPI
rank count according to the node size, for example:

```
export OMP_NUM_THREADS=32
bin/gendie_preprocessor case -f EXAMPLE -n 32 --skip-modules
```

You can also run the example from the source tree after building, but the driver
expects the required executables next to `preprocessor.py`. Installing first is
the simpler and less error-prone path for new users.

The driver prints stage markers as it runs:

```
[0] start
[1] preprocessing configuration
[2] mindist
[3] coarse mesher
[4] mesh filter
[5] monitor function
[6] mesh refinement
[7] fine mesh filter
```

## Parallelization
Parts of the application are parallelized using either OpenMP or MPI. Set a
sane value for `OMP_NUM_THREADS` before running `gendie_preprocessor` and pass
the desired number of MPI ranks via `-n`.

The OpenMP and MPI phases run sequentially, not at the same time. This means a
dedicated single-node run can use the full node width for both
`OMP_NUM_THREADS` and `-n`. See
[`docs/SingleNodeClusterGuide.md`](docs/SingleNodeClusterGuide.md) for details.



Example setup.e3d for Box-Mesh
```
[E3DGeometryData]
[E3DGeometryData/Preprocessing]
HexMesher=Box
sEl_x = 1.0
sEl_y = 0.6
sEl_z = 1.0
geometryStart = -330.00,-104.50,0.00
geometryLength = +660.0,209.0,520.0
```

Example setup.e3d for HollowCylinder-Mesh
```
[E3DGeometryData]
[E3DGeometryData/Preprocessing]
HexMesher=HollowCylinder
sEl_Tangential = 1.25
sEl_Radial = 0.66
sEl_Axial = 1.0
BarrelDiameter = 1000.0
InnerDiameter = 440.0
BarrelLength = 625.0
AxialStartPosition = -54.0
```

Example setup.e3d for FullCylinder-Mesh
```
[E3DGeometryData]
[E3DGeometryData/Preprocessing]
HexMesher=FullCylinder
FullCylinderPeriodicity = 4   # optional, defaults to 4 when omitted
sEl_Tangential = 1.25
sEl_Radial = 0.80
sEl_Axial = 1.25
BarrelDiameter = 1000.0
BarrelLength = 600.0
AxialStartPosition = -54.0
```
