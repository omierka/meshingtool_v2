#pragma once

#include <cgal_types.hpp>
#include <meshhexer/types.hpp>

namespace MeshHexer
{
  /**
   * \brief Compute the diameter of the maximal inscribed sphere at each face of the mesh.
   *
   * The maximum inscribed sphere of a face is the largest possible sphere
   * that touches the centroid of the face and any other point of the mesh, without intersecting the mesh.
   *
   * The diameter of each sphere is made available in a "f:MIS_diameter" mesh property of type double.
   * The id of the _other_ primitive that is touched by the sphere is made available in a "f:MIS_id" mesh property.
   *
   * See the following for details on the algorithm:
   * Shrinking sphere: A parallel algorithm for computing the thickness of 3D objects
   * Masatomo Inui, Nobuyuki Umezu, Ryohei Shimane
   * COMPUTER-AIDED DESIGN & APPLICATIONS, 2016, VOL. 13, NO. 2, 199–207
   * http://dx.doi.org/10.1080/16864360.2015.1084186
   */
  void maximal_inscribed_spheres(Mesh& mesh, const AABBTree& aabb_tree);

  /**
   * \brief Calculates topological distances on the mesh
   *
   * \param mesh Mesh to calculate distances on
   * \param property Name of a mesh property of with type FaceIndex
   *
   * Calculates the smallest distance along the edges of the mesh between any face f
   * and the corresponding face property[f].
   */
  void topological_distances(Mesh& mesh, const std::string& property, double max_distance = 0.0);
  void topological_distances(Mesh& mesh, const std::string& targets_property, const std::string& max_distance_property);
  void update_validity_from_neighbor_diameters(Mesh& mesh);
  void update_validity_from_neighbor_normals(Mesh& mesh);
  void update_validity_from_small_angles(Mesh& mesh, double min_angle_degrees = 4.0);

  /**
   * \brief Computes vertex normals and makes them available as v:normals
   */
  void compute_vertex_normals(Mesh& mesh);
  void compute_max_dihedral_angle(Mesh& mesh);
  void compute_curvature(Mesh& mesh);

  /**
   * \brief Returns true if face normals calculated via a cross product point towards the "unbounded" side of the mesh
   *
   * We write "unbounded" because we do not assume that the mesh is watertight. There thus might not be a "bounded"
   * side.angle This function will work either way.
   */
  bool do_normals_point_outside(const Mesh& mesh);

  bool is_wound_consistently(const Mesh& mesh);

  /**
   * \brief Returns largest length of AABB surrounding the mesh
   */
  double mesh_size(const Mesh& mesh);
  BoundingBox bounding_box(const Mesh& mesh);

  void score_gaps(Mesh& mesh);

  std::vector<Gap> gaps(Mesh& mesh);

  /// Determine a best-effort guess at the minimal gap of the given mesh
  Gap min_gap(Mesh& mesh);

  /// Shoot inward normals per face and store nearest hit distance / target face
  void normal_distances(Mesh& mesh, const AABBTree& aabb_tree);

  /// Combine MIS and normal distances into a single monitor field
  void monitor_distances(Mesh& mesh);
  void invalidate_monitors_below(Mesh& mesh, double min_gap_diameter);
  struct MonitorHistogramBin
  {
    double lower;
    double upper;
    double area;
  };
  struct MonitorHistogramComponentStats
  {
    std::size_t collections = 0;
    double connected_area = 0.0;
    double isolated_area = 0.0;
    double max_component_area = 0.0;
  };
  std::vector<MonitorHistogramBin> monitor_histogram(Mesh& mesh, double min_gap_diameter);
  std::vector<MonitorHistogramComponentStats> monitor_histogram_component_counts(
    Mesh& mesh,
    const std::vector<MonitorHistogramBin>& bins,
    double min_gap_diameter,
    double min_component_area_fraction);
  double adjusted_min_gap_from_histogram(const std::vector<MonitorHistogramBin>& bins, double min_gap_diameter);
  std::pair<std::size_t, std::size_t> monitor_histogram_min_max_indices(const std::vector<MonitorHistogramBin>& bins);
  /// Faces whose centroid normals fall below the requested magnitude
  std::vector<FaceIndex> faces_with_short_vertex_normals(Mesh& mesh, double min_length);
  /// Faces whose geometric area falls below a relative threshold
  std::vector<FaceIndex> faces_with_small_area(const Mesh& mesh, double relative_threshold);

  std::vector<std::pair<Point2D, double>> z_depths(Mesh& mesh, AABBTree& aabb_tree);
} // namespace MeshHexer
