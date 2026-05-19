# Single-Node Cluster Guide

This note is intended for cluster deployment of the New Gen Meshing Tool when the full workflow is always executed on exactly one compute node.

## Practical Command

For the bundled example case, run from the repository root:

```bash
python3 preprocessor.py case -f EXAMPLE -n <mpi_ranks>
```

After installation, run from the install root:

```bash
gendie_preprocessor case -f EXAMPLE -n <mpi_ranks>
```

`-n` controls the MPI rank count used by the MPI-based stages.

## Workflow Sequence

The driver executes these components in order:

1. `meshhexer-cli report --configure-case-for-preprocessing <case>`
2. `meshhexer-cli min-gap <case>/surface.off`
3. `generate_hollow_cylinder_mesh -i <case>/setup.e3d -o <case>/Coarse_meshDir/Mesh.tri -s <coarse_size>`
4. `mpirun -np <mpi_ranks> meshcleaner ...`
5. `mpirun -np <mpi_ranks> hex_VS_triangulation_intersection ...`
6. `meshref -f <case> -d <span>`
7. `mpirun -np <mpi_ranks> meshcleaner ...`

If the monitor post-check detects that the first volumetric monitor bin is effectively empty, the driver applies one correction step and repeats the meshing part of the workflow. In that case, stages 3 to 7 run a second time.

## Parallel Model By Component

| Step | Component | Parallel model | Notes |
|---|---|---|---|
| 1 | `meshhexer-cli report` | OpenMP | Built with `MESHHEXER_HAVE_OMP=ON`. |
| 2 | `meshhexer-cli min-gap` | OpenMP | Same executable as step 1. |
| 3 | `generate_hollow_cylinder_mesh` | Serial | Python helper, no MPI or OpenMP. |
| 4 | `meshcleaner` | MPI | Launched through `mpirun`. Source uses `MPI_Init`. |
| 5 | `hex_VS_triangulation_intersection` | MPI | Launched through `mpirun`. Source uses `MPI_Init`. |
| 6 | `meshref` | Serial | No MPI and no OpenMP in the current implementation. |
| 7 | `meshcleaner` | MPI | Same executable as step 4. |

There is no stage in the current `preprocessor.py` workflow that uses MPI and OpenMP simultaneously inside the same executable. The workflow mixes OpenMP and MPI across different sequential stages.

## Single-Node Scheduling Rule

Run the job inside an allocation that contains exactly one compute node. All MPI ranks must stay on that one node.

Because the stages are sequential, there is no runtime overlap between:

- the OpenMP stages (`meshhexer-cli`)
- the MPI stages (`meshcleaner`, `hex_VS_triangulation_intersection`)

That means the node can be sized once and then reused by each stage in turn.

## Recommended Runtime Configuration

Use these settings as the default single-node configuration:

```bash
export OMP_NUM_THREADS=<omp_threads>
export OMP_PROC_BIND=close
export OMP_PLACES=cores
python3 preprocessor.py case -f EXAMPLE -n <mpi_ranks>
```

Recommended choices:

- choose `<mpi_ranks>` as the number of physical CPU cores you want to give to the MPI stages on the node
- choose `<omp_threads>` as the number of CPU cores you want `meshhexer-cli` to use in the OpenMP stages
- keep both values less than or equal to the number of physical cores on the node

Since the OpenMP and MPI phases do not overlap, it is acceptable on a dedicated node to use the full node width for both, for example:

```bash
export OMP_NUM_THREADS=64
python3 preprocessor.py case -f EXAMPLE -n 64
```

This does not create 64 x 64 concurrent workers, because the OpenMP stages finish before the MPI stages start.

## Conservative Starting Point

If Dirk wants a low-risk first configuration for a fresh cluster environment:

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
python3 preprocessor.py case -f EXAMPLE -n 32
```

Then tune based on which phase dominates wall-clock time:

- if `meshhexer-cli` is slow, increase `OMP_NUM_THREADS`
- if `meshcleaner` or `hex_VS_triangulation_intersection` is slow, increase `-n`

## Important MPI Detail

`meshcleaner` is MPI-based and its own documentation requires at least two MPI ranks. In practice, use:

```bash
-n 2
```

as the minimum valid cluster setting, and typically much higher on a full compute node.

## Example Slurm Job

```bash
#!/bin/bash
#SBATCH --job-name=ngmt-example
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

module purge
module load cmake/3.28.3 gcc/latest-v13 openmpi/4.1.6 python/3.13.5

export OMP_NUM_THREADS=32
export OMP_PROC_BIND=close
export OMP_PLACES=cores

cd /path/to/newgenmeshingtool
python3 preprocessor.py case -f EXAMPLE -n 32
```

This layout keeps the run on one node, gives the MPI stages 32 ranks, and gives the OpenMP stages 32 threads.

## Operational Summary

- `meshhexer-cli` uses OpenMP.
- `meshcleaner` uses MPI.
- `hex_VS_triangulation_intersection` uses MPI.
- `meshref` is serial.
- `generate_hollow_cylinder_mesh` is serial.
- The workflow is mixed-model overall, but each stage is either OpenMP, MPI, or serial.
- For cluster use, reserve one node and tune `OMP_NUM_THREADS` and `-n` independently.
