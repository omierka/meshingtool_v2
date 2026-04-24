# Mesh Analysis README

This document describes the round-geometry analysis currently available in `MeshHexer`. The feature is intended for `FullCylinder` and `HollowCylinder` surface triangulations and is used through `meshhexer-cli report --configure-case-for-preprocessing`.

## 1. Purpose

The analysis computes geometric quantities needed later in the preprocessing and meshing workflow:

- xy offset of the cylindrical axis with respect to `(0, 0)`
- outer diameter of the geometry
- inner diameter of the geometry
- classification as `FullCylinder` or `HollowCylinder`
- common inlet/outlet extrusion length
- physical axial extent `z_min_phys`, `z_max_phys`

The implementation assumes:

- the main cylinder axis is parallel to the global `z` axis
- tilt is not allowed
- inflow normals point into the closed triangulation
- the geometry is a regular round geometry for which global plane sweeps do not hit unrelated features in misleading ways

## 2. Entry Point

Run:

```bash
meshhexer-cli report --configure-case-for-preprocessing <case-folder>
```

The command expects a `setup.e3d` file in the same folder as the mesh:

```text
<folder>/surface.off
<folder>/setup.e3d
```

From `setup.e3d` the analysis currently reads, for each inflow:

- `center`
- `normal`

These values are read from sections of the form:

```text
[E3DProcessParameters/Inflow_<N>]
center = x,y,z
normal = nx,ny,nz
```

## 3. Axis Offset

The axis offset is determined from a horizontal slice near the top of the geometry:

```text
z_slice = z_max - 1e-4 * (z_max - z_min)
```

Steps:

1. Intersect the triangulation with the plane `z = z_slice`.
2. Keep the largest closed intersection loop in xy projection.
3. Fit a circle to that loop.
4. Use the fitted center as the cylinder-axis estimate `(x_c, y_c)`.

The geometry is reported as aligned to the origin if:

```text
abs(x_c) <= 1e-4 and abs(y_c) <= 1e-4
```

This tolerance is read from `preprocessor.cfg` via:

```ini
[MeshAnalysis]
axis_alignment_epsilon = 1.0e-4
```

## 4. Outer Diameter

The outer diameter is computed from the full triangulation, not from a single slice.

Steps:

1. Use the fitted axis center `(x_c, y_c)`.
2. For every vertex of the triangulation compute the radial distance

```text
r = sqrt((x - x_c)^2 + (y - y_c)^2)
```

3. Take the global maximum radial distance as the raw outer radius.
4. Apply inflow-based restrictions.

### Inflow restriction for the outer radius

If an inflow normal points toward the axis, that inflow limits the usable outer radius.

For an inflow center `c` and inflow normal `n`, the radial direction `e_r` is taken from the fitted axis to the inflow center. If the radial projection satisfies

```text
dot(n_xy_normalized, e_r) < 0
```

then the inflow center radius is treated as an upper bound for the outer radius.

The final outer radius is therefore:

```text
outer_radius = min(global_outer_radius, all_outer_inflow_limits)
```

and the reported outer diameter is:

```text
outer_diameter = 2 * outer_radius
```

## 5. Inner Diameter

The inner diameter is also computed from the full triangulation, using the complete surface as a global constraint.

The intended geometric meaning is:

- the largest infinite cylinder
- aligned with the global `z` axis
- centered at the fitted axis `(x_c, y_c)`
- that can be placed inside the triangulation without intersecting the surface

Implementation:

1. Project every triangle of the surface triangulation to the xy plane.
2. Compute the distance from the fitted axis point `(x_c, y_c)` to each projected triangle.
3. Take the minimum such distance over the complete triangulation.

That minimum radial distance is the inner radius:

```text
inner_radius = min distance from (x_c, y_c) to any projected surface triangle
```

and:

```text
inner_diameter = 2 * inner_radius
```

This definition makes the inner diameter a global feasibility criterion for an endless z-aligned cylinder. If the axis is blocked anywhere by the triangulation, the inner diameter becomes zero.

## 6. Full vs Hollow Classification

After computing outer and inner diameters, the geometry is classified by the ratio:

```text
inner_diameter / outer_diameter
```

Current rule:

- if the ratio is `< 0.10`, classify as `FullCylinder`
- otherwise classify as `HollowCylinder`

The threshold is currently hard-coded in the analysis config and is intended to become configurable later.

## 7. Extrusion Length

The feature also estimates the common artificial inlet/outlet extrusion length from the inflows.

For each inflow:

1. Start at the inflow center.
2. Move conceptually opposite to the inflow normal.
3. Determine how far the triangulation extends in that direction.

Implementation detail:

- for each inflow, all mesh vertices are projected onto the inflow normal
- the most negative projection relative to the inflow center determines the extrusion length sample

So for one inflow:

```text
L_ext = max(0, -min(dot(p - c_inflow, n_inflow)))
```

For multiple inflows:

- one sample is computed per inflow
- the median sample is reported as `extrusion_length`
- all samples are checked against a relative consistency tolerance of `5 %`

The report therefore contains:

- `extrusion_length`
- `extrusion_length_consistent`
- `extrusion_length_samples`

## 8. Physical z Extent

The triangulation may be axially longer than the physical CFD domain because of artificial extrusions.

The analysis reports:

- `z_min_phys`
- `z_max_phys`

using:

```text
z_max_phys = z_max_tri - extrusion_length
```

The `z_min` side is shortened only if there is an inflow whose normal qualifies as an axial `z-` limiter.

Current rule:

- normalize the inflow normal
- if `n_z > cos(5 deg)`, then that inflow is treated as limiting the `z_min` side

Then:

```text
z_min_phys = z_min_tri + extrusion_length
```

otherwise:

```text
z_min_phys = z_min_tri
```

## 9. Reported Quantities

The `report --configure-case-for-preprocessing` path currently prints:

- `Axis center (x, y)`
- `Axis aligned to origin`
- `Top slice z`
- `Outer diameter`
- `Inner diameter`
- `Inner/outer ratio`
- `Classification`
- `Extrusion length`
- `Extrusion length consistent`
- `z-min limited by axial inflow`
- `z_min_phys`
- `z_max_phys`
- `Extrusion samples`

## 10. Example

```bash
meshhexer-cli report --configure-case-for-preprocessing 15992
```

This uses:

- `15992/surface.off`
- `15992/setup.e3d`

and prints the normal mesh report together with an additional `Round Geometry Analysis` section.

## 11. Current Limitations

- The feature currently lives behind `report --configure-case-for-preprocessing`; there is no separate machine-readable command yet.
- The implementation currently parses only `center` and `normal` from `setup.e3d`.
- The analysis is written for round geometries only. Box analysis is not part of this feature yet.
- The global extrusion-length logic intentionally treats any global intersection as relevant; this is acceptable only under the assumption of regular round geometries.
- The axis-alignment epsilon is exposed through `preprocessor.cfg` under `[MeshAnalysis]`.
- Thresholds such as the `10 %` full/hollow criterion and the `5 deg` axial inflow criterion are not yet exposed through `preprocessor.cfg`.

## 12. Recommended Usage

Use this feature:

- for `FullCylinder` and `HollowCylinder` cases only
- when a matching `setup.e3d` with inflow definitions is present next to the surface mesh

Do not use it as a generic geometry classifier for arbitrary triangulations.
