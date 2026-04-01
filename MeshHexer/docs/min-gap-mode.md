# Min-Gap Mode Workflow

This document summarises the processing pipeline behind the `min-gap` mode of `meshhexer-cli` together with the artefacts it produces. It is intended for users who need to understand which geometric quantities are generated, how they are combined, and how the defaults can be overridden as the workflow evolves.

## 1. Entry Point

```
meshhexer-cli --checkpoint-path <out/MINGAP.vtu> min-gap <mesh/surface.off>
```

The CLI expects a triangulated surface mesh (`.off`, `.ply`, `.obj`, …). The optional `--checkpoint-path` flag instructs MeshHexer to persist every intermediate mesh property (`MIS_diameter`, `normaldistance`, `Monitor`, `Validity`, …) to a VTU file that downstream tools can post-process.

The command prints **three** whitespace-separated values:

1. `mindist` – the adjusted minimal gap (see Section 5)
2. `span` – an integer in `1..3` denoting how wide the dominant monitor distribution is in logarithmic bins
3. `coarseMeshSize` – `mindist * coarse_mesh_scaling * 3^span`, with `coarse_mesh_scaling` defaulting to 1.3

All three scalars are consumed by `runner.sh`. The driver now forwards the emitted coarse mesh size directly to the structured mesher, so changing the scaling factor is centralised inside the min-gap workflow.

## 2. Geometry Preparation

Before any gap is scored the code ensures that the following fields are available on the mesh:

| Field                     | Purpose                                                                |
|---------------------------|------------------------------------------------------------------------|
| `v:normals`               | Vertex normals used by several validity tests                          |
| `f:MIS_diameter`, `f:MIS_id` | Maximal inscribed sphere diameter and limiting face id              |
| `f:normaldistance`, `f:normaldistance_target` | Inward ray-cast distance to the closest opposite triangle |
| `f:Validity`              | Integer flag; values 1–5 mark triangles rejected by various heuristics |
| `f:topological_distance`  | Path length along the surface between a face and its MIS partner       |
| `f:gap_score`             | Confidence that a given MIS really captures a physical gap             |

During this stage MeshHexer also applies three validity heuristics:

1. **Neighbour diameters** – if all neighbours are >2× larger the face is flagged (2)
2. **Neighbour normals** – >60° difference flags the face (3)
3. **Small angles** – any interior angle <4° flags the face (4)

Faces touching their own MIS partner receive flag (1). All flagged faces still carry actual values in the VTU file but are excluded from the min-gap decision.

## 3. Monitor Field

The monitor value of face _f_ is defined as `max(MIS_diameter[f], normaldistance[f])`. This scalar approximates the “practical” gap likely to be hit by a CFD mesh. The full map is written into the VTU and later filtered as described in Section 5.

## 4. Histogram

After the raw min-gap (based on MIS diameters only) has been found, MeshHexer builds an area-weighted histogram of all monitor values that survive the validity checks.

- The first bin spans `[mindist, 3¹·mindist]`, the second `[3¹·mindist, 3²·mindist]`, etc.  
- All bins exist in the exported `size_distribution_histogram.txt`, even when cumulative fractions are tiny.
- Each row contains:
  - `Mark`: `"MIN"`, `"MAX"`, `"MIN/MAX"`, or empty
  - `Bin Start`, `Bin End`: lower/upper thresholds
  - `Area`: sum of triangle areas inside the bin
  - `Percentage`: share of total weighted area

### Selecting MIN/MAX bins

1. Starting at the first bin, accumulate percentages until the sum exceeds 0.1 %. The bin where this happens is tagged `MIN`.
2. From that bin upward, keep accumulating until reaching 80 % coverage or until three additional bins have been included. The last contributing bin is tagged `MAX`. If the mesh only has a single populated bin, it receives the `MIN/MAX` mark.

The integer `span = max_index − min_index` (clamped to `[1, 3]`) corresponds to the width column returned by the CLI.

## 5. Adjusted Min-Gap

After emitting the histogram MeshHexer adjusts the scalar `mindist` by setting it to the lower edge of the `MIN` bin. All monitor values strictly lower than that threshold are overwritten with `-1` and assigned validity flag 5, so post-processing tools can easily mask them out.

This updated `mindist` (together with the histogram span and the derived coarse mesh size) is what the CLI prints and what the runner uses to create subsequent coarse meshes.

## 6. Produced Artefacts

Running the `min-gap` pipeline with `--checkpoint-path ${FOLDER}/MINGAP.vtu` produces:

| File                                | Description                                                         |
|-------------------------------------|---------------------------------------------------------------------|
| `${FOLDER}/MINGAP.vtu`              | Surface mesh annotated with all intermediate properties             |
| `${FOLDER}/size_distribution_histogram.txt` | Area-weighted monitor histogram with MIN/MAX markers          |
| CLI stdout (`mindist span coarseMeshSize`) | Scalars ingested by `runner.sh`                               |

Downstream stages use these artefacts for mesh cleaning, monitor-field creation, and refinement, with `CoarseMeshSize = mindist * coarse_mesh_scaling * 3^span` (default scaling 1.3) providing a conservative starting resolution.

## 7. Adjusting Defaults

Most thresholds mentioned above are hard-coded today (0.1 %, 66 %, 4°…). When changes are needed:

1. Update the corresponding helper in `src/properties.cpp` (e.g. `update_validity_from_small_angles` or `monitor_histogram_min_max_indices`).
2. Rebuild `meshhexer-cli`.
3. The CLI/runner wiring will automatically pick up the new behaviour, as only the values flowing through the existing pipeline change.

Please coordinate changes to the console output (`mindist span coarseMeshSize`) with users of `runner.sh` to keep the coarse-mesh derivation consistent.

## 8. Prechecking and Repairing Meshes

The min-gap pipeline terminates early when it encounters invalid triangles (zero-length centroid normals, degenerate or near-zero-area faces). Before running Section 1 make sure the input mesh passes:

```
# Exit status 0 means “safe for min-gap”
meshhexer-cli precheck <mesh/surface.off>
```

`precheck` lists up to `--max-report` problematic faces and exits with code 2 if any are found, so runners can gate their workflows accordingly.

If `precheck` reports degenerate/near-zero-area faces you can repair the surface directly:

```
# Writes mesh_repaired.off unless --output is given
meshhexer-cli repair <mesh/surface.off>
meshhexer-cli precheck <mesh_repaired.off>   # should now pass
```

`repair` removes every flagged face, deletes isolated vertices, and triangulates the remaining surface so the cleaned mesh can be fed back into the min-gap workflow.
