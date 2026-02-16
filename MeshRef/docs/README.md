# Hex Mesh Refinement – Developer Notes

## Repository Layout

```
Makefile                Build rules for the `meshref` executable.
src/
  var.f90               Global template catalogue, mesh data types, helpers.
  inout.f90             File I/O for templates/target meshes and VTU writers.
  def.f90               Refinement logic, randomization (pure mesh processing).
  cleanup.f90           Patch deduplication (intra- and inter-element) + refined mesh assembly.
main.f90              CLI front-end (`-f`, `-r`, `-t` options).
VERTEX/, PATCHES/, ...  Sample template meshes provided with the project.
generate_tri_files.py   Utility for creating TRI meshes from canonical patterns.
doc/README.md           This document.
decode_templates.py     Helper tool for mapping 8-bit patterns to template IDs.
```

The code assumes a Fortran compiler such as `gfortran`. Build with `make` (override `FC` if needed) which produces the `meshref` binary in the project root. Run `./meshref -f <working-folder> [ -r <0-100> ] [ -t <threshold> ]`. The working folder must contain `Coarse_meshDir/Mesh.tri`, `area.txt`, and optionally `setup.e3d`. The executable writes `target.vtu`, `RefinedCleanMesh_lvl0*.vtu`, and `RefinedCleanMesh_lvl1*.vtu` inside that folder (raw, `_clean`, and `_refined_clean` variants) and, after the last cleanup pass, emits a TRI file at `<folder>/meshDir_BU/Merged_Mesh.tri`. If `setup.e3d` specifies `HexMesher=HollowCylinder` under `[E3DGeometryData/Preprocessing]`, the cylindrical coordinate transformation path is activated automatically and a message is printed to the terminal.

The `-t/--threshold` option lets you override the monitor cut-off (`Monitor_threshold` in `var.f90`) that separates the randomly marked elements. When omitted, the default of `1.5` is used.

## Module Overview

### `var_mod` (`src/var.f90`)

This module contains:

* Template catalogue: `mesh_file_count`, the list of canonical template filenames, and `templates(8,mesh_file_count)` holding the 8‑bit (0/1) patterns. `template_is_final(:)` marks whether a template is terminal (`FIN/`) or requires further consideration (`INT/`).
* Global mesh storage:
* `type(mesh_type)` holds TRI data (counts, coordinates, connectivity, original `KNPR` flags, the per-element span lists `kelementspan(:)`, and the per-element `monitor(:)` array used for CLI-driven random refinement). Instances include `meshes(:)` (templates) and the hierarchy `hex_mesh(0:depth)` that stores every refinement level (`hex_mesh(0)` is the original target, `hex_mesh(1)` the first cleaned refinement, `hex_mesh(2)` the inter-element-cleaned mesh ready for the next stage, etc.). `initialize_mesh_levels(depth)`/`release_mesh_levels()` manage this hierarchy, and `bind_refined_mesh(level)` updates the module pointers `target_mesh`/`refined_clean_mesh`. The module also surfaces tunables such as `Monitor_threshold` and `scaling_factor_TRI_output` (defaults to `0.1`) used when exporting TRI backups.
  * `type(element_patch_type)` captures a refined template embedded into target coordinates (local/global coords, connectivity, per-vertex `KNPR`). `element_patch_group` stores the final patches per target element (`patchlist(:)`, `count`), and `element_patches(:)` is an array of these groups. `clean_element_patches(:)` mirrors this structure and, after cleanup, holds deduplicated node/element data per root element (single patch per group). `refined_clean_mesh` (also `type(mesh_type)`) stores the globally deduplicated coordinate/connectivity arrays produced by the inter-patch cleanup stage.
* Helper routines:
  * `clear_mesh`, `clear_patch`, `release_element_patches`.
  * Pattern utilities: extract template code strings, convert to logical masks, rotate via `rotate_patch`, compare via `determine_template`.
* `match_template_pattern` uses filename digits and the canonical table to locate the appropriate template ID.
* `clear_patch_group` / `append_patch_to_group` manage the per-element patch lists.

### `inout_mod` (`src/inout.f90`)

Responsibilities:

* `load_all_meshes` walks the canonical template list, reads each TRI file, resolves its 8‑bit code, logs the mapping, and stores the geometry/connectivity into `meshes`. It also marks `template_is_final(idx)` based on the FIN/INT folder.
* `load_target_mesh` reads `Coarse_meshDir/Mesh.tri` from the CLI-specified folder. After loading, it inspects `setup.e3d` for `HexMesher=HollowCylinder` and toggles cylindrical template insertion accordingly.
* `read_mesh_file` parses the custom TRI format: skips headers, reads counts, coordinate block (`DCORVG`), connectivity block (`KVERT`), and nodal property block (`KNPR`). All arrays are allocated to exact sizes.
* `write_patch_group_vtu` takes any array of `element_patch_group` structures (either the raw patches or the cleaned ones), concatenates their nodes/cells into a VTU file, and emits both the per-vertex `KNPR` data and a per-cell `monitor` integer array. Helpers such as `build_clean_output_path`, `build_refined_clean_output_path`, and `build_level_output_path` centralize the filename manipulation needed for each refinement level. `write_refined_clean_tri` mirrors the structure of the input TRI file (header, `DCORVG`, `KVERT`, `KNPR`) and writes the global refined mesh to `meshDir_BU/Merged_Mesh.tri`, scaling coordinates by `scaling_factor_TRI_output` (default `0.1`) and setting all exported `KNPR` values to zero.

### `def_mod` (`src/def.f90`)

High-level controls:

* `initialize_mesh_database`/`finalize_mesh_database` prepare and release the template catalogue and optional target data.
* `random_element_marking(percent)` (CLI option `-r`) populates the per-element `monitor(:)` array with random 0/1/2 tags (selected elements randomly choose between 1 and 2, unselected stay 0). A module-level flag `reproducibility` (currently default `.true.`) keeps the pseudo-random streams deterministic across runs; flip it to `.false.` if you prefer unique randomness per execution.
* `vertice_marking(min_marker)` converts the `monitor(:)` field into vertex `KNPR` flags by marking every vertex of elements whose monitor value is at least `min_marker`.
* `enforce_refinement_levels(mesh, level)` walks the element adjacency lists (`kelementspan`) and ensures that neighbors of any element marked at `level` are promoted to at least `level - 1`, smoothing the refinement field prior to recursion.

Mesh refinement workflow:

1. `mesh_refinement(mesh)`
   * Validates that the target mesh and template database are loaded.
   * Allocates one `element_patch_group` per target element (cleared each run).
   * Derives the per-element neighbor lists (`kelementspan`) for downstream processing.
   * For each target element, calls the recursive worker `refine_element`, which keeps refining non-final templates until either `%is_final` is true or the depth cap `max_refinement_depth` (currently `1`) is reached.
   * Tracks template usage/unknown patterns and, once recursion finishes, prints the summary. VTU output is intentionally handled by the caller now. The routine accepts any `mesh_type` (e.g., `target_mesh` or `refined_clean_mesh`) so future refinement passes can re-enter on previously cleaned meshes.

2. `refine_element`
   * Builds the logical pattern from the current element’s `KNPR`, rotates it to canonical orientation, and determines the template ID.
   * Uses `instantiate_patch` to embed the template in the parent coordinate system, calling `fill_up_element` to generate global coordinates.
   * If the patch is non-final and the depth limit hasn’t been met, every sub-element undergoes the same process. Otherwise the patch is stored permanently via `store_final_patch`.

3. `fill_up_element`
   * Computes the trilinear mapping between the canonical reference cube and the target element using the provided 8 corner coordinates.
   * Fills `patch%global_coor(:, :)` with the transformed positions, ensuring the first eight entries align exactly with the original element vertices (preserving ordering for downstream refinement and VTU export).

### `cleanup_mod` (`src/cleanup.f90`)

* `clear_intra_patches` consumes the raw `element_patches` array, builds a per-root merged patch in `clean_element_patches`, and prints temporary diagnostics showing vertex counts before/after deduplication. The tolerance used for deduplication is measured upfront (half the shortest hexa edge from the supplied patches, or `1e-10` fallback); the main loop also calls this estimator before invoking the cleanup so the user sees it in the logs.
* `clear_inter_patches([span_mesh])` walks `clean_element_patches`, uses the supplied mesh’s (default `target_mesh`) `kelementspan` incidence to compare only neighboring root elements, and constructs `refined_clean_mesh` with a globally deduplicated `coor/kvert` pair (one shared vertex set for all refined cells); VTU export is handled separately via `write_refined_clean_vtu`.
* All coordinates within a root element are deduplicated with a small tolerance and the `kvert` connectivity is remapped to the new unique node indices. `clean_element_patches(i)%patchlist` therefore contains at most one aggregated patch per root element, and `monitor` values are preserved for every refined sub-element.
* Before any cleanup runs, the module measures the shortest edge found across all refined hexahedra (from the raw `element_patches`) and sets the merge tolerance to half that length, falling back to `1e-10` if no edges are available.

### `main.f90`

`meshref_main` is the program entry point. It parses the command line:

* `-f, --folder <dir>` (mandatory) – working folder; expects `Coarse_meshDir/Mesh.tri`, `area.txt`, and optionally `setup.e3d`.
* `-r, --random-refinement <0-100>` (optional) – randomize the target `KNPR` before refinement (overrides `area.txt`).
* `-t, --threshold <value>` (optional) – override the monitor threshold used when interpreting the monitor distribution.
* `-h, --help` – display usage.

Execution order:

1. Parse arguments.
2. Initialize templates (`initialize_mesh_database`) and load the target mesh.
3. Optionally call `random_element_marking` followed by `vertice_marking(1)` (controlled by `-r`). If `--folder` supplies an `area.txt`, those values override the random marks. The resulting monitor field is then “smoothed” via `enforce_refinement_levels`, which ensures that an element marked for the highest refinement level is surrounded only by elements of the same level or exactly one level lower.
4. For each refinement level `L` (currently hard-coded to two levels: `L = 0` and `L = 1`), run `mesh_refinement(hex_mesh(L))`, write the raw VTU (`*_lvlL.vtu`), run intra-patch cleanup and emit `*_lvlL_clean.vtu`, run inter-patch cleanup using the corresponding span data, build the next level’s span for `hex_mesh(L+1)`, and write the inter-element VTU (`*_lvlL_refined_clean.vtu`).
5. Finalize and free resources.

## Typical Workflow

1. `make` – builds `meshref`.
2. `./meshref -f <working-folder> [ -r <0-100> ] [ -t <threshold> ]` – loads canonical templates (logged to stdout), reads `Coarse_meshDir/Mesh.tri`, writes `target.vtu` (the unrefined mesh carrying the current `monitor`/`KNPR` state), recursively refines each element (respecting the depth cap), and writes `RefinedCleanMesh_lvl0*.vtu`/`RefinedCleanMesh_lvl1*.vtu` (raw, `_clean`, `_refined_clean`). After the final inter-element cleanup the assembled mesh is also exported as `<folder>/meshDir_BU/Merged_Mesh.tri` using the original header metadata but updated `NEL/NVT` counts. The run also assembles `hex_mesh(1)` and `hex_mesh(2)` in memory for future inter-element connectivity processing.
3. The program writes `*_lvlL.vtu` containing the raw element patches for each level. It then runs the cleanup stage, emitting a second file whose name inserts `_clean` (e.g., `*_lvlL_clean.vtu`). This cleaned VTU reuses the aggregated per-root patches with deduplicated coordinates/connectivity, making intra-patch connectivity explicit. Finally, `clear_inter_patches` merges coincident nodes across adjacent root elements, populates the next entry in `hex_mesh`, and `write_refined_clean_vtu` emits `*_lvlL_refined_clean.vtu`.
4. Open either VTU in ParaView for inspection. Both carry the propagated `KNPR` flags on `PointData`; the cleaned file has fewer duplicate vertices per root element.

## Extending the Code

* New canonical templates: add the TRI file to `mesh_files`, add its 8‑bit pattern to `templates`, and ensure `template_is_final` is set correctly (the loader also tags FIN/INT automatically based on folder names).
* Additional diagnostics or refinement strategies hook naturally into `mesh_refinement` by inspecting the per-element patterns before constructing patches.
* The VTU writer can be extended with additional `PointData` or `CellData` arrays by mimicking the existing `KNPR` output block.

## Utilities

`decode_templates.py` accepts one or more 8-bit combinations (e.g. `11101000` or `[10111110]`) and reports the canonical template ID/code after applying the same rotation logic used by `var_mod`. Run `./decode_templates.py <pattern> [...]`, provide a filename that lists codes (e.g. `./decode_templates.py list`), or pass `--file list.txt` for batch decoding. Each output line shows the original digits, the canonical code, and whether the match was direct or required rotation.

`generate_tri_files.py` (moved from `TESTCASES/`) generates synthetic TRI meshes that follow the expected VERTEX directory layout. Use it to create new template or test meshes from prescribed coordinate/connectivity inputs without editing TRI files by hand.

For further tweaks (e.g., deeper recursion, custom refinement rules), `refine_element`, `instantiate_patch`, and `fill_up_element` are the primary entry points to intercept and transform template instances before export.
