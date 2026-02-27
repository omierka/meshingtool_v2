codex resume 019c8a82-f198-7ea0-90ee-dd64e4712593

MeshDeform
==========

MeshDeform applies a deformation monitor to the filtered hexahedral mesh that comes out of the standard workflow (`Filtered.tri`).  It loads the hexahedral mesh and the surrounding `surface.off` triangulation referenced in `setup.e3d`, assigns parametrization constraints to every boundary vertex, and then performs one or more umbrella-style deformation steps that push the mesh toward the signed-distance monitor while respecting the boundary parameterizations.

### High-level workflow

1. **Input parsing (`meshdeform_main`)**
   * reads CLI arguments (`--folder` is mandatory)
   * loads `setup.e3d` to get box/cylinder parameters and process inflows
   * loads the `surface.off` triangulation and the `Filtered.tri` mesh
2. **Preprocessing (`def_mod`)**
   * builds topological boundary masks (`identify_boundary_nodes`)
   * tags vertices with parametrization bits (`assign_vertex_constraints`)
   * mirrors mesh data and the CGAL surface mesh to every MPI rank (`broadcast_mesh`, `broadcast_mesh_coordinates`)
3. **Distance estimation**
   * each rank computes signed distances for its subset of vertices using `compute_signed_distances_parallel`, which calls the CGAL AABB queries and gathers the result with `MPI_Allgatherv`
4. **Deformation step**
   * rank 0 runs `apply_edge_deformation`, which:
     * computes volume-weighted umbrella averages along the hexahedral edges
     * enforces boundary parametrizations (plane/cylinder/corner constraints)
     * blends the new position with damping (currently 0.2)
   * updated coordinates are broadcast for the next distance sweep
5. **Output**
   * rank 0 writes `MeshDeformResult.vtu` (includes the signed-distance field) and `Mesh_deformed.tri`

### Key routines

| Routine | Description |
| --- | --- |
| `compute_signed_distances_parallel` | MPI-parallel CGAL distance evaluation (per-vertex signed distance to `surface.off`) |
| `apply_edge_deformation` | Performs one umbrella deformation step, respecting boundary parametrizations |
| `write_deformed_vtu` / `write_deformed_tri` | Writes the updated mesh (VTU with distance data, TRI with new coordinates) |
| `assign_vertex_constraints` / `enforce_parametrization_constraints` | Classify boundary vertices and restrict motion to planes/cylinders |

### Configuration knobs

`preprocessor.cfg` now has a `[3DMeshDeform]` section:

```
[3DMeshDeform]
deformation_steps = 1
characteristic_scale = 100.0
```

* `deformation_steps` controls how many umbrella iterations run per invocation.
* `characteristic_scale` scales the signed-distance weighting relative to the characteristic mesh size.

### Usage

```
mpirun -np <N> meshdeform --folder <case-folder>
```

* `--folder` must point to a case directory containing `Filtered.tri`, `surface.off`, and `setup.e3d`.
* `--input` can override the TRI file (defaults to `<folder>/Filtered.tri`).
* `--verbose` prints extra configuration and monitor details on rank 0.

The program prints a progress banner such as

```
Performing deformation steps : [1][2]...[N]
Deformation time (s):   3.42
```

All deformation work happens on rank 0; helper ranks are used solely for the CGAL distance evaluation, so `-np 1` runs remain valid if MPI is unavailable.***
