#include <cgal_types.hpp>
#include <limits>
#include <macros.hpp>
#include <properties.hpp>
#include <meshhexer/types.hpp>

#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/heap/priority_queue.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <omp.h>

namespace MeshHexer
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  namespace
  {
    Vector3D surface_normal(Mesh& mesh, FaceIndex f, Point3D point)
    {
      auto maybe_normals_map = mesh.property_map<VertexIndex, Vector3D>("v:normals");
      XASSERTM(maybe_normals_map.has_value(), "No vertex normals!");

      const Mesh::Property_map<VertexIndex, Vector3D>& vertex_normals = maybe_normals_map.value();

      // Get barycentric coordinates of point
      auto location = CGAL::Polygon_mesh_processing::locate_in_face(point, f, mesh);

      std::array<Vector3D, 3> normals;

      // Get vertex normals
      int i(0);
      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(f)))
      {
        normals.at(i++) = vertex_normals[v];
      }

      if(i == 3)
      {
        return location.second[0] * normals[0] + location.second[1] * normals[1] + location.second[2] * normals[2];
      }
      return {0.0, 0.0, 0.0};
    }

    // Set up priority queue for choosing next vertex to explore
    struct FrontierEntry
    {
      VertexIndex idx;
      double priority;

      FrontierEntry() : idx(0), priority(0)
      {
      }

      FrontierEntry(VertexIndex v, double p) : idx(v), priority(p)
      {
      }
    };

    // The priority_queue returns the last element of the defined ordering first.
    // We want the last element to be the one with the lowest estimated distance.
    // Thus we sort in descending order.
    struct Comparator
    {
      bool operator()(FrontierEntry a, FrontierEntry b) const
      {
        return a.priority > b.priority;
      }
    };

    double topological_distance(
      FaceIndex a,
      FaceIndex b,
      const Mesh& mesh,
      const std::vector<double>& edge_lengths,
      double max_distance = 0.0)
    {
      // Allocate frontier
      boost::heap::priority_queue<FrontierEntry, boost::heap::compare<Comparator>> frontier;

      // Allocate distances
      boost::unordered_flat_map<VertexIndex, double> distances;

      // Determine end points
      std::vector<VertexIndex> goal_vertices;
      std::vector<Point3D> goal_points;
      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(b)))
      {
        goal_vertices.push_back(v);
        goal_points.push_back(mesh.point(v));
      }

      // Minimum direct distance between any start and end point.
      // If this distance is larger then the maximum search distance
      // we do not have to start any path-finding, because the
      // points can not be connected within the max distance
      double minimum_distance(std::numeric_limits<double>::max());

      // Iterate over starting points
      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(a)))
      {
        for(const Point3D& goal : goal_points)
        {
          minimum_distance = std::min(minimum_distance, (mesh.point(v) - goal).squared_length());
        }

        // Add start point to frontier
        frontier.emplace(v, 0);

        // Start point has best known distance 0.0
        distances[v] = 0.0;
      }

      if(max_distance > 0.0 && std::sqrt(minimum_distance) > max_distance)
      {
        return -1.0;
      }

      // Set up heuristic
      // Heuristic is shortest direct distance to any of the goal vertices
      const Point3D heuristic_target_point = goal_points.front();
      const auto heuristic = [&](VertexIndex v)
      {
        return std::sqrt((mesh.point(v) - heuristic_target_point).squared_length());
      };

      while(!frontier.empty())
      {
        // Remove next best element from the frontier
        VertexIndex current = frontier.top().idx;
        frontier.pop();

        double current_distance = distances[current];

        if(std::find(goal_vertices.begin(), goal_vertices.end(), current) != goal_vertices.end())
        {
          // We found one of the target vertices.
          return current_distance;
        }

        for(HalfedgeIndex hedge : mesh.halfedges_around_target(mesh.halfedge(current)))
        {
          VertexIndex neighbor = mesh.source(hedge);

          double edge_length = edge_lengths[mesh.edge(hedge)];
          double new_cost = current_distance + edge_length;

          if(max_distance != 0.0 && new_cost > max_distance)
          {
            continue;
          }

          auto [it, inserted] = distances.try_emplace(neighbor, std::numeric_limits<double>::max());

          if(new_cost < it->second)
          {
            it->second = new_cost;
            double priority = new_cost + heuristic(neighbor);
            frontier.emplace(neighbor, priority);
          }
        }
      }

      return -1.0;
    }
  } // namespace

  BoundingBox bounding_box(const Mesh& mesh)
  {
    const Point3D& p = mesh.point(VertexIndex(0));
    Point min = {p.x(), p.y(), p.z()};
    Point max = {p.x(), p.y(), p.z()};
    for(const auto& vertex : mesh.points())
    {
      min.x = std::min(min.x, vertex.x());
      min.y = std::min(min.y, vertex.y());
      min.z = std::min(min.z, vertex.z());

      max.x = std::max(max.x, vertex.x());
      max.y = std::max(max.y, vertex.y());
      max.z = std::max(max.z, vertex.z());
    }

    return BoundingBox{min, max};
  }

  double mesh_size(const Mesh& mesh)
  {
    BoundingBox bb = bounding_box(mesh);
    return std::max({bb.max.x - bb.min.x, bb.max.y - bb.min.y, bb.max.z - bb.min.z});
  }

  bool is_wound_consistently(const Mesh& mesh)
  {
    for(HalfedgeIndex idx : mesh.halfedges())
    {
      VertexIndex source = mesh.source(idx);
      VertexIndex target = mesh.target(idx);

      HalfedgeIndex opposite = mesh.opposite(idx);

      if(opposite != Mesh::null_halfedge())
      {
        // If the mesh is wound consistently, then opposite halfedges need
        // to visit the same vertices in the opposite order.
        // Otherwise the two faces are wound differently.
        if(source != mesh.target(opposite) || target != mesh.source(opposite))
        {
          return false;
        }
      }
    }
    return true;
  }

  // Compute the maximum dihedral angle for each face in radians
  void compute_max_dihedral_angle(Mesh& mesh)
  {
    Mesh::Property_map<FaceIndex, double> dihedral_angles =
      mesh.add_property_map<FaceIndex, double>("f:dihedral_angle", 0).first;

    MESHHEXER_PRAGMA_OMP(parallel for)
    for(FaceIndex current : mesh.faces())
    {
      double max_angle = 0;
      Vector3D normal = PMP::compute_face_normal(current, mesh);

      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(current)))
      {
        for(FaceIndex neighbor : mesh.faces_around_target(mesh.halfedge(v)))
        {
          if(neighbor > mesh.num_faces())
          {
            // NOTE(mmuegge): The inner loop produces invalid indices for the 10062 mesh
            continue;
          }
          Vector3D n = PMP::compute_face_normal(neighbor, mesh);

          double cosine = CGAL::scalar_product(normal, n);
          cosine = std::max(std::min(cosine, 1.0), 0.0);

          max_angle = std::max(max_angle, std::acos(cosine));
        }
      }

      dihedral_angles.begin();

      dihedral_angles[current] = max_angle;
    }
  }

  void compute_vertex_normals(Mesh& mesh)
  {
    Mesh::Property_map<VertexIndex, Vector3D> normals =
      mesh.add_property_map<VertexIndex, Vector3D>("v:normals", Vector3D(0.0, 0.0, 0.0)).first;

    CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, normals);
  }

  void compute_curvature(Mesh& mesh)
  {
    Mesh::Property_map<VertexIndex, double> mean_curvature =
      mesh.add_property_map<VertexIndex, double>("v:mean_curvature", 0).first;

    Mesh::Property_map<VertexIndex, double> gaussian_curvature =
      mesh.add_property_map<VertexIndex, double>("v:gaussian_curvature", 0).first;

    Mesh::Property_map<VertexIndex, PMP::Principal_curvatures_and_directions<Kernel>> principal_curvatures =
      mesh
        .add_property_map<VertexIndex, PMP::Principal_curvatures_and_directions<Kernel>>(
          "v:principal_curvatures",
          PMP::Principal_curvatures_and_directions<Kernel>(
            Real(0),
            Real(0),
            Vector3D(0.0, 0.0, 0.0),
            Vector3D(0.0, 0.0, 0.0)))
        .first;

    PMP::interpolated_corrected_curvatures(
      mesh,
      CGAL::parameters::vertex_mean_curvature_map(mean_curvature)
        .vertex_Gaussian_curvature(gaussian_curvature)
        .vertex_principal_curvatures_and_directions_map(principal_curvatures));
  }

  void maximal_inscribed_spheres(Mesh& mesh, const AABBTree& aabb_tree)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    if(mesh.num_faces() == 0)
    {
      // No work to be done
      return;
    }

    // Create mesh property for storing thicknesses
    Mesh::Property_map<FaceIndex, double> diameter_property =
      mesh.add_property_map<FaceIndex, double>("f:MIS_diameter", 0).first;

    Mesh::Property_map<FaceIndex, std::uint32_t> id_property =
      mesh.add_property_map<FaceIndex, std::uint32_t>("f:MIS_id", 0).first;

    Mesh::Property_map<FaceIndex, double> similarity =
      mesh.add_property_map<FaceIndex, double>("f:similarity_of_normals", 0).first;

    Mesh::Property_map<FaceIndex, int> iters = mesh.add_property_map<FaceIndex, int>("f:MIS_iters", 0).first;

    // Determine bounding box of mesh
    const Real max_radius = std::min({
                              aabb_tree.bbox().xmax() - aabb_tree.bbox().xmin(),
                              aabb_tree.bbox().ymax() - aabb_tree.bbox().ymin(),
                              aabb_tree.bbox().zmax() - aabb_tree.bbox().zmin(),
                            }) /
                            Real(2.0);

    const Real normal_direction = PMP::is_outward_oriented(mesh) ? Real(-1.0) : Real(1.0);

    MESHHEXER_PRAGMA_OMP(parallel for)
    for(FaceIndex face_index : mesh.faces())
    {
      // Determine centroid of face
      Point3D centroid(0.0, 0.0, 0.0);
      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(face_index)))
      {
        centroid += Real(1.0 / 3.0) * Vector3D(Point3D(CGAL::Origin()), mesh.point(v));
      }

      // Determine average edge length around the face.
      // We later use this to judge if a new closest point is
      // actually a closer point or just a slightly closer
      // vertex of the discretization of the same inscribed sphere

      Real local_max_edge_length(0);
      for(HalfedgeIndex e : mesh.halfedges_around_face(mesh.halfedge(face_index)))
      {
        local_max_edge_length = std::max(local_max_edge_length, PMP::edge_length(e, mesh));
      }

      // Determine initial sphere radius
      // We need an initial radius that is large enough for our initial sphere
      // to intersect at least one object other than the current face.
      // We cast a ray from the centroid towards the face normal.
      // If we hit something, we have an initial guess for the radius.
      // Otherwise the current face is pointing towards a hole in the surface mesh.
      // In that case we use the smallest length of the axis aligned bounding box as the initial diameter.
      // This is the largest value any MIS can be.

      Real initial_radius(max_radius);
      FaceIndex initial_id(0);

      // Determine (inward) normal
      // Optimization(mmuegge): We know the barycentric coordinates of the centroid already.
      // No reason for surface_normal to recalculate them
      // Either offer a version of surface_normal that accepts a point in barycentric coordinates
      // or a centroid_normal function.
      Vector3D normal(surface_normal(mesh, face_index, centroid));
      normal /= std::sqrt(normal.squared_length());
      Vector3D inward_normal = normal_direction * normal;

      // Cast ray
      Ray3D ray(centroid, inward_normal);
      auto skip = [=](FaceIndex idx) { return idx == face_index; };
      std::optional<RayIntersection> intersection = aabb_tree.first_intersection(ray, skip);

      if(intersection && std::holds_alternative<Point3D>(intersection.value().first))
      {
        const Point3D& p = std::get<Point3D>(intersection.value().first);
        Real distance = CGAL::approximate_sqrt(Vector3D(centroid, p).squared_length());

        // NOTE(mmuegge): In meshes with self-intersections it can occur that a raycast hits
        // a triangle distinct from the current face at distance 0.
        // In that case we want to keep the initial radius based on the bounding box.
        if(distance > 0.0)
        {
          initial_radius = std::min(initial_radius, distance / 2.0);
          initial_id = intersection.value().second;
        }
      }

      if(intersection && std::holds_alternative<Segment3D>(intersection.value().first))
      {
        const Segment3D& segment = std::get<Segment3D>(intersection.value().first);
        Real distance = CGAL::approximate_sqrt(Vector3D(centroid, segment.source()).squared_length());

        if(distance > 0.0)
        {
          initial_radius = std::min(initial_radius, distance / 2.0);
          initial_id = intersection.value().second;
        }
      }

      Real prev_radius = initial_radius + 1;
      Real radius = initial_radius;
      FaceIndex id = initial_id;
      Point3D closest_point;

      // Shrink sphere until change in radius becomes too small
      while(prev_radius - radius > 1e-6)
      {
        Point3D center = centroid + radius * inward_normal;
        auto closest_point_data = aabb_tree.closest_point_and_primitive(center);

        if(closest_point_data.second == face_index)
        {
          // Disregard closest points on self
          break;
        }

        closest_point = closest_point_data.first;

        Vector3D to_closest = Vector3D(centroid, closest_point);
        if(CGAL::scalar_product(inward_normal, to_closest) <= Real(0))
        {
          // Next closest point is behind the current face.
          // Disregard it
          break;
        }

        Real distance_to_closest = CGAL::approximate_sqrt(Vector3D(center, closest_point).squared_length());

        Real discretization_error(0);
        // NOTE(mmuegge): Must be local_max_edge_length < 2 * radius to be a valid triangle
        // Tests show 1.5 works better though
        if(local_max_edge_length < 1.5 * radius)
        {
          // Largest distance between sphere and an edge of maximum length, assuming the vertices of that edge lie on
          // the sphere. Any closest points within this distance belong to the same radius, they are just slightly
          // offset due to the surface discretization.
          discretization_error =
            radius - CGAL::approximate_sqrt((radius * radius) - std::pow(local_max_edge_length / 2.0, 2));
        }
        else
        {
          discretization_error = 0.02 * radius;
        }

        if(distance_to_closest >= radius - 2 * discretization_error)
        {
          // Closest point is not in sphere. Current sphere is correct
          break;
        }

        // We now need to find the radius r' of the new sphere that touches both the
        // centroid of our face and the closest point we just found.
        // Let the centroid be p, the closest point p' and let c' be the center of the new sphere.
        // p, p', and c' prime form a isosceles triangle with two sides of length r', a side of length d.
        //
        //     r'
        // p'------c'
        // \       |
        //  \      |
        //   \     |
        //  d \    | r'
        //     \   |
        //      \  |
        //       \a|
        //        p
        //
        // The angle a is the angle between the normal of the face and the vector pp'.
        // We can then determine r' via the law of sines as r' = d * sin(a) * sin(180 - 2a)

        double d = CGAL::to_double(CGAL::approximate_sqrt((closest_point - centroid).squared_length()));
        double alpha =
          CGAL::to_double(CGAL::approximate_angle(inward_normal, Vector3D(centroid, closest_point))) / 180.0 * M_PI;

        Real new_radius = NAN;
        if(alpha != 0)
        {
          new_radius = Real(d * std::sin(alpha) / std::sin(M_PI - (2.0 * alpha)));
        }
        else
        {
          new_radius = d;
        }

        prev_radius = radius;
        radius = std::min(new_radius, radius);
        id = closest_point_data.second;

        iters[face_index]++;
      }

      // Record thickness
      diameter_property[face_index] = CGAL::to_double(2.0 * radius);
      id_property[face_index] = id;

      // Record abs(cos(theta))
      // similarity[face_index] = std::abs(CGAL::to_double(CGAL::scalar_product(inward_normal, surface_normal(mesh, id,
      // closest_point))));
      similarity[face_index] = CGAL::scalar_product(inward_normal, -PMP::compute_face_normal(id, mesh));
    }
  }

  void topological_distances(Mesh& mesh, const std::string& targets_property, double max_distance)
  {
    // Get targets
    auto maybe_target_map = mesh.property_map<FaceIndex, std::uint32_t>(targets_property);
    if(!maybe_target_map.has_value())
    {
      return;
    }
    Mesh::Property_map<FaceIndex, std::uint32_t> targets = maybe_target_map.value();

    // Get output property
    Mesh::Property_map<FaceIndex, double> topo_distance =
      mesh.add_property_map<FaceIndex, double>("f:topological_distance", 0).first;

    // Pre-calculate edge lengths
    std::vector<double> edge_lengths;
    for(EdgeIndex idx : mesh.edges())
    {
      Point3D a = mesh.point(mesh.target(mesh.halfedge(idx)));
      Point3D b = mesh.point(mesh.source(mesh.halfedge(idx)));
      edge_lengths.push_back(std::sqrt((a - b).squared_length()));
    }

    MESHHEXER_PRAGMA_OMP(parallel for schedule(dynamic))
    for(FaceIndex f : mesh.faces())
    {
      topo_distance[f] = topological_distance(f, FaceIndex(targets[f]), mesh, edge_lengths, max_distance);
    }
  }

  void topological_distances(Mesh& mesh, const std::string& targets_property, const std::string& max_distance_property)
  {
    // Get targets
    auto maybe_target_map = mesh.property_map<FaceIndex, std::uint32_t>(targets_property);
    if(!maybe_target_map.has_value())
    {
      return;
    }
    Mesh::Property_map<FaceIndex, std::uint32_t> targets = maybe_target_map.value();

    // Get maximum search distances
    auto maybe_max_distances_map = mesh.property_map<FaceIndex, double>(max_distance_property);
    if(!maybe_max_distances_map.has_value())
    {
      return;
    }
    Mesh::Property_map<FaceIndex, double> max_distances = maybe_max_distances_map.value();

    // Get output property
    Mesh::Property_map<FaceIndex, double> topo_distance =
      mesh.add_property_map<FaceIndex, double>("f:topological_distance", 0).first;

    // Pre-calculate edge lengths
    std::vector<double> edge_lengths;
    for(EdgeIndex idx : mesh.edges())
    {
      Point3D a = mesh.point(mesh.target(mesh.halfedge(idx)));
      Point3D b = mesh.point(mesh.source(mesh.halfedge(idx)));
      edge_lengths.push_back(std::sqrt((a - b).squared_length()));
    }

    MESHHEXER_PRAGMA_OMP(parallel for schedule(dynamic))
    for(FaceIndex f : mesh.faces())
    {
      topo_distance[f] = topological_distance(f, FaceIndex(targets[f]), mesh, edge_lengths, max_distances[f]);
    }
  }

  void score_gaps(Mesh& mesh)
  {
    MeshHexer::Mesh::Property_map<MeshHexer::FaceIndex, double> diameters =
      mesh.property_map<MeshHexer::FaceIndex, double>("f:MIS_diameter").value();

    MeshHexer::Mesh::Property_map<MeshHexer::FaceIndex, double> topo_dists =
      mesh.property_map<MeshHexer::FaceIndex, double>("f:topological_distance").value();

    MeshHexer::Mesh::Property_map<MeshHexer::FaceIndex, double> similarity_of_normals =
      mesh.property_map<MeshHexer::FaceIndex, double>("f:similarity_of_normals").value();

    MeshHexer::Mesh::Property_map<MeshHexer::FaceIndex, double> dihedral_angles =
      mesh.property_map<MeshHexer::FaceIndex, double>("f:dihedral_angle").value();

    MeshHexer::Mesh::Property_map<MeshHexer::FaceIndex, std::uint32_t> ids =
      mesh.property_map<MeshHexer::FaceIndex, std::uint32_t>("f:MIS_id").value();

    MeshHexer::Mesh::Property_map<MeshHexer::FaceIndex, double> max_search_distances =
      mesh.property_map<MeshHexer::FaceIndex, double>("f:max_search_distance").value();

    auto topo_dist_score = [&](FaceIndex f)
    {
      double dist = topo_dists[f];

      // If the topological_distance is zero, this is not a valid gap
      if(dist == 0.0)
      {
        return 0.0;
      }

      // If the faces are either unconnected or farther apart than the search radius,
      // then this is an ideal gap
      if(dist == -1.0)
      {
        return 1.0;
      }

      // Divide by maximum search range for this face
      double rel = topo_dists[f] / max_search_distances[f];

      // Relative distance to opposite side of sphere
      double opposite = ((M_PI / 2.0) * diameters[f]) / max_search_distances[f];

      if(opposite < 1.0 && rel < opposite)
      {
        // Opposite side of sphere is reachable within max_search_distance and limiting face is
        // closer than opposite side.
        // Empirically scoring quadratically works well. We thus want the score to be
        // a * rel^2 with a * opposite^2 = 1 => a = 1 / opposite^2
        return (1.0 / (opposite * opposite)) * rel * rel;
      }

      // Above a relative distance of opposite there must be some amount of space between the MIS and the mesh.
      // All these gaps are ok.
      return 1.0;
    };

    auto similarity_of_normals_score = [&](FaceIndex f)
    {
      if(-similarity_of_normals[f] < 0.0)
      {
        return 0.0;
      }

      return std::pow(-similarity_of_normals[f], 1.0 / 3.0);
    };

    auto dihedral_angle_score = [&](FaceIndex f)
    {
      if(dihedral_angles[f] < 0.3)
      {
        return 1.0;
      }
      return std::sqrt(1.0 - ((dihedral_angles[f] + 0.3) / M_PI));
    };

    auto max_edge_length = [&](FaceIndex f)
    {
      double result = 0.0;
      for(HalfedgeIndex idx : mesh.halfedges_around_face(mesh.halfedge(f)))
      {
        result = std::max(result, PMP::edge_length(idx, mesh));
      }
      return result;
    };

    std::vector<std::pair<MeshHexer::FaceIndex, MeshHexer::FaceIndex>> self_intersections;
    CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(self_intersections));

    Mesh::Property_map<FaceIndex, double> gap_score = mesh.add_property_map<FaceIndex, double>("f:gap_score", 0).first;

    double max_diameter = mesh_size(mesh);

    for(FaceIndex f : mesh.faces())
    {
      auto limiting = FaceIndex(ids[f]);

      double edge_ratio = max_edge_length(f) / max_edge_length(limiting);
      if(edge_ratio < 1.0)
      {
        edge_ratio = 1.0 / edge_ratio;
      }

      if(
        PMP::face_aspect_ratio(f, mesh) < 5.0 && PMP::face_aspect_ratio(limiting, mesh) < 5.0 &&
        std::find(self_intersections.begin(), self_intersections.end(), std::make_pair(f, limiting)) ==
          self_intersections.end() &&
        std::find(self_intersections.begin(), self_intersections.end(), std::make_pair(limiting, f)) ==
          self_intersections.end() &&
        diameters[f] / max_diameter > 1e-4 && edge_ratio < 5.0)
      {
        if(topo_dists[f] == -1.0 && topo_dists[FaceIndex(ids[f])] == -1.0)
        {
          // Perfect gap if both triangles are unconnected.
          // NOTE(mmuegge): There are situations where a large triangle
          // hits a very small triangle across a corner of the mesh.
          // In that case the small triangle might not find its limiting triangle
          // within the search radius for topological distances, but the large
          // triangle might, because its vertices are further from its centroid,
          // giving an initial distance boost.
          // We thus check that neither triangle found its limiting triangle,
          // to determine if we should score the gap more exactly.
          // Note that the above condition does not imply f = ids[f].
          // TODO(mmuegge): Can we fix this by considering the distances
          // from the centroid to a starting vertex in the topological distance search?
          gap_score[f] = 1.0;
        }
        else
        {
          gap_score[f] = (1.0 / 3.0) * (topo_dist_score(f) + similarity_of_normals_score(f) + dihedral_angle_score(f));
        }
      }
      else
      {
        gap_score[f] = 0;
      }
    }
  }

  std::vector<Gap> gaps(Mesh& mesh)
  {
    std::vector<Gap> result;

    Mesh::Property_map<FaceIndex, double> gap_scores =
      mesh.property_map<FaceIndex, double>("f:gap_score").value();

    Mesh::Property_map<FaceIndex, std::uint32_t> ids =
      mesh.property_map<FaceIndex, std::uint32_t>("f:MIS_id").value();

    Mesh::Property_map<FaceIndex, double> gaps =
      mesh.property_map<FaceIndex, double>("f:MIS_diameter").value();

    for(const FaceIndex f : mesh.faces())
    {
      Point3D centroid(0.0, 0.0, 0.0);
      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(f)))
      {
        centroid += Real(1.0 / 3.0) * Vector3D(Point3D(CGAL::Origin()), mesh.point(v));
      }

      result.emplace_back(static_cast<std::uint32_t>(f), ids[f], gaps[f], gap_scores[f]);
    }

    return result;
  }

  Gap min_gap(Mesh& mesh)
  {
    MeshHexer::Mesh::Property_map<MeshHexer::FaceIndex, double> diameters =
      mesh.property_map<MeshHexer::FaceIndex, double>("f:MIS_diameter").value();

    MeshHexer::Mesh::Property_map<MeshHexer::FaceIndex, double> scores =
      mesh.property_map<MeshHexer::FaceIndex, double>("f:gap_score").value();

    MeshHexer::Mesh::Property_map<MeshHexer::FaceIndex, std::uint32_t> ids =
      mesh.property_map<MeshHexer::FaceIndex, std::uint32_t>("f:MIS_id").value();

    double threshold = 0.95;

    std::vector<FaceIndex> candidates;
    std::copy_if(
      mesh.faces_begin(),
      mesh.faces_end(),
      std::back_inserter(candidates),
      [&](FaceIndex f) { return scores[f] > threshold; });

    if(candidates.empty())
    {
      // No gaps with score above 0.95. Take 10th percentile instead.

      std::vector<double> percentiles(scores.begin(), scores.end());
      // Sort scores in descending order
      std::sort(percentiles.begin(), percentiles.end(), [](double a, double b) { return a > b; });

      std::size_t percentile_size = std::ceil(Real(percentiles.size()) * Real(0.1));
      double cutoff_score = *(percentiles.begin() + percentile_size);

      std::copy_if(
        mesh.faces_begin(),
        mesh.faces_end(),
        std::back_inserter(candidates),
        [&](FaceIndex f) { return scores[f] > cutoff_score; });
    }

    FaceIndex gap_face(0);
    double minimum = std::numeric_limits<double>::max();

    for(FaceIndex f : candidates)
    {
      if(diameters[f] < minimum)
      {
        gap_face = f;
        minimum = diameters[f];
      }
    }

    return {gap_face, ids[gap_face], diameters[gap_face], scores[gap_face]};
  }

  std::vector<std::pair<Point2D, double>> z_depths(Mesh& mesh, AABBTree& aabb_tree)
  {
    const Real start_z = aabb_tree.bbox().zmin() - Real(1);
    const Vector3D ray_dir(0.0, 0.0, 1.0);

    std::vector<std::pair<Point2D, double>> result;

    // Shoot a ray from start_z in direction -z through all vertices of the mesh

    std::vector<RayIntersection> intersections;

    for(const FaceIndex face : mesh.faces())
    {
      if(CGAL::scalar_product(PMP::compute_face_normal(face, mesh), ray_dir) >= -1e-4)
      {
        // Face is orthogonal to ray direction or points away. Skip it.
        continue;
      }

      Point3D centroid(0.0, 0.0, 0.0);
      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(face)))
      {
        centroid += Real(1.0 / 3.0) * Vector3D(Point3D(CGAL::Origin()), mesh.point(v));
      }

      intersections.clear();

      const Point3D ray_origin(centroid.x(), centroid.y(), start_z);
      const Ray3D ray(ray_origin, ray_dir);

      aabb_tree.all_intersections(ray, std::back_inserter(intersections));

      // Walk through intersections
      // Every second invertval is inside the mesh and counts as depth

      auto it = intersections.begin();

      double depth = 0;
      while(it + 1 < intersections.end())
      {
        RayIntersection& entry = *it;
        RayIntersection& exit = *(it + 1);

        it += 2;

        if(!std::holds_alternative<Point3D>(entry.first) || !std::holds_alternative<Point3D>(exit.first))
        {
          // Ray passes along face, rather than through it.
          // Just ignore this face for now
          break;
        }

        const Point3D& entry_point = std::get<Point3D>(entry.first);
        const Point3D& exit_point = std::get<Point3D>(exit.first);

        depth += exit_point.z() - entry_point.z();
      }

      if(depth > aabb_tree.bbox().zmax() - aabb_tree.bbox().zmin() + 1)
      {
        break;
      }

      result.emplace_back(Point2D(ray_origin.x(), ray_origin.y()), depth);
    }

    return result;
  }
} // namespace MeshHexer
