# New Gen Meshing Tool

## Dependencies
Required:
- C++ and Fortran compilers
- CMake (3.28.0 or greater)
- meshref (from FeatFlower)

Optional, will be downloaded if missing:
- CGAL (6.0.0 or greater)
- Boost (1.88 or greater)

## Building

Extract the NewGenMeshingTool archive.
Create and enter build directory:

```
mkdir path/to/build/directoy && cd path/to/build/directory
```

Call CMake to set up the build system:

```
cmake path/to/root/of/archive
```

Then build the project via:
```
cmake --build .
```

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
to install the tool to a custom location. Then add the custom location to your ```PATH``` via
```
export PATH=/path/to/install/location:$PATH
```

The tool requires both itself and ```meshref``` to be in the PATH.

## Usage
Run
```
newgenmeshingtool --folder /path/to/input/directory
```
to generate a coarse mesh.

The input directory given as an argument to ```newgenmeshingtool``` needs to contain the following files:
- ```surface.off```, the surface mesh to generate a coarse mesh for
- ```setup.e3d```, preprocessing and simulation settings
- ```param_meshref.cfg```,

The directory you are calling ```newgenmeshingtool``` from needs to contain the following files:
- a ```PATCHES``` directory, containing the necessary patches for ```meshref```
- a directory ```start``` containing a ```sampleRigidBody.xml```

## Parallelization
Parts of the application are parallelized using either OMP or MPI. Set a sane value for ```OMP_NUM_THREADS``` before running ```newgenmeshingtool``` and ensure it runs on a machine that can allocate 64 mpi processes.



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
