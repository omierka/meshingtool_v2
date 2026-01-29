#include <cgal_types.hpp>
#include <meshhexer/types.hpp>
#include <meshing.hpp>
#include <properties.hpp>

#include <cstdint>
#include <limits>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/enum.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Octree.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/property_map.h>

namespace MeshHexer
{
  namespace PS = CGAL::Polyline_simplification_2;

  using Octree = CGAL::Octree<
    Kernel,
    std::vector<std::pair<Point3D, Point3D>>,
    CGAL::First_of_pair_property_map<std::pair<Point3D, Point3D>>>;

  namespace
  {
    enum class WindingOrder : std::uint8_t
    {
      Clockwise,
      CounterClockwise,
    };

    /**
     * \brief Compute left normal of line segment from a to b. Result has unit length.
     */
    Vector2D left_normal(const Point2D& a, const Point2D& b)
    {
      Vector2D delta = b - a;
      Vector2D normal = Vector2D(-delta.y(), delta.x());
      return normal / CGAL::approximate_sqrt(normal.squared_length());
    }

    /**
     * \brief Compute right normal of line segment from a to b. Result has unit length.
     */
    Vector2D right_normal(const Point2D& a, const Point2D& b)
    {
      return -left_normal(a, b);
    }

    Vector2D outside_normal(const Point2D& a, const Point2D& b, WindingOrder order)
    {
      return order == WindingOrder::CounterClockwise ? right_normal(a, b) : left_normal(a, b);
    }

    double angle(const Vector2D& a, const Vector2D& b)
    {
      return CGAL::to_double(CGAL::approximate_angle(Vector3D(a.x(), a.y(), 0), Vector3D(b.x(), b.y(), 0)));
    }

    /**
     * \brief Returns true if all vertices of \c b are contained in \c a.
     */
    bool contains(const Polygon2D& a, const Polygon2D& b)
    {
      return std::all_of(
        b.begin(),
        b.end(),
        [&a](const Point2D& p) { return a.bounded_side(p) != CGAL::ON_UNBOUNDED_SIDE; });
    }

    std::vector<PolygonWithHoles2D> make_polygon(Polylines2D& polylines)
    {
      // Create polygons from the polylines.
      // This gives us inside-out tests without
      // without having to think about the orientation
      // of the polyline
      std::vector<Polygon2D> polygons;

      for(auto& polyline : polylines)
      {
        Polygon2D poly(polyline.begin(), polyline.end() - 1);
        if(poly.is_clockwise_oriented())
        {
          polygons.emplace_back(std::reverse_iterator(polyline.end() - 1), std::reverse_iterator(polyline.begin()));
        }
        else
        {
          polygons.emplace_back(std::move(poly));
        }
      }

      // Simplyify polygons
      for(Polygon2D& poly : polygons)
      {
        poly = PS::simplify(poly, PS::Squared_distance_cost(), PS::Stop_above_cost_threshold(1e-4));
      }

      // Sort polygons in descending order by area, i.e. largest polygon comes first.
      std::sort(
        polygons.begin(),
        polygons.end(),
        [](const Polygon2D& a, const Polygon2D& b) { return a.area() > b.area(); });

      // Create inclusion tree for polygons.
      const std::size_t num_polygons = polygons.size();
      std::vector<int> depths(num_polygons, 0);
      std::vector<std::size_t> parents(num_polygons, -1);

      std::size_t current_idx = 0;
      for(const Polygon2D& polygon : polygons)
      {
        std::size_t parent_candidate = ~std::size_t(0);
        int parent_depth = -1;
        // Find the deepest, i.e. smallest, already handled polygon
        // that still contains the current polygon
        for(std::size_t j(0); j < current_idx; j++)
        {
          if(contains(polygons[j], polygon) && depths[j] > parent_depth)
          {
            parent_depth = depths[j];
            parent_candidate = j;
          }
        }

        if(parent_candidate != ~std::size_t(0))
        {
          // Polygon is contained in another polygon. Set parent
          parents[current_idx] = parent_candidate;
          depths[current_idx] = parent_depth + 1;
        }
        else
        {
          // Polygon is not contained in another larger polygon. Add it as new root
          parents[current_idx] = current_idx;
          depths[current_idx] = 0;
        }
        current_idx++;
      }

      std::vector<std::size_t> root_mapping(num_polygons, -1);
      std::vector<PolygonWithHoles2D> result;

      std::size_t handled_polygons = 0;

      // Collect all roots
      for(std::size_t i(0); i < num_polygons; i++)
      {
        if(depths[i] == 0)
        {
          root_mapping[i] = result.size();
          result.emplace_back(polygons[i]);
          handled_polygons++;
        }
      }

      // Collect holes at depth 1
      for(std::size_t i(0); i < num_polygons; i++)
      {
        if(depths[i] == 1)
        {
          std::size_t parent = parents[i];
          result[root_mapping[parent]].add_hole(polygons[i]);
          handled_polygons++;
        }
      }

      std::size_t unhandled_polygons = num_polygons - handled_polygons;
      if(unhandled_polygons > 0)
      {
        std::cerr << "Warning: " << unhandled_polygons << " unhandled polygons in MeshHexer::make_polygon!";
      }

      return result;
    }

    bool merge_polygons_with_holes(const PolygonWithHoles2D& a, const PolygonWithHoles2D& b, PolygonWithHoles2D& res)
    {
      if(CGAL::join(a.outer_boundary(), b.outer_boundary(), res))
      {
        std::vector<PolygonWithHoles2D> intersections;

        for(const auto& hole_a : a.holes())
        {
          for(const auto& hole_b : b.holes())
          {
            CGAL::intersection(hole_a, hole_b, std::back_inserter(intersections));
          }
        }

        std::vector<Polygon2D> new_holes;

        new_holes.reserve(intersections.size());
        for(auto& intersection : intersections)
        {
          new_holes.push_back(intersection.outer_boundary());
        }

        res.holes().clear();
        res.holes().insert(res.holes().end(), new_holes.begin(), new_holes.end());
        return true;
      }
      else
      {
        return false;
      }
    }

    std::vector<PolygonWithHoles2D> merge(std::vector<PolygonWithHoles2D>& cross_sections)
    {
      std::vector<PolygonWithHoles2D> result;
      while(!cross_sections.empty())
      {
        PolygonWithHoles2D merged = cross_sections.back();
        cross_sections.pop_back();

        bool progress = false;
        do // NOLINT
        {
          progress = false;
          for(auto it = cross_sections.rbegin(); it != cross_sections.rend();)
          {
            if(it->outer_boundary().is_clockwise_oriented())
            {
              it = decltype(it){cross_sections.erase(std::next(it).base())};
              continue;
            }

            PolygonWithHoles2D new_shadow;
            if(merge_polygons_with_holes(merged, *it, new_shadow))
            {
              it = decltype(it){cross_sections.erase(std::next(it).base())};

              merged = new_shadow;
              progress = true;
              continue;
            }
            it++;
          }
        } while(progress);

        merged = PS::simplify(merged, PS::Squared_distance_cost(), PS::Stop_above_cost_threshold(1e-8));
        result.push_back(merged);
      }

      return result;
    }

    template<typename OutputIterator>
    void
    find_cross_section(CGAL::Polygon_mesh_slicer<Mesh, Kernel>& slicer, const CuttingPlane& plane, OutputIterator out)
    {
      Polylines3D polylines;
      Polylines2D projected_polylines;

      // Find polylines
      polylines.clear();
      projected_polylines.clear();
      slicer(plane.plane, std::back_inserter(polylines));

      if(polylines.empty())
      {
        std::cerr << "Empty intersection.\n";
        return;
      }

      // The found polylines are 3d objects. We are actually only interested in the projection onto the current cutting
      // plane.
      for(const auto& polyline : polylines)
      {
        Polyline2D projected;
        for(const Point3D& point : polyline)
        {
          projected.push_back(plane.project(point));
        }
        projected_polylines.push_back(projected);
      }

      std::vector<PolygonWithHoles2D> new_cross_sections = make_polygon(projected_polylines);

      for(PolygonWithHoles2D& poly : new_cross_sections)
      {
        *out = poly;
      }
    }

    template<typename OutputIterator>
    void find_cross_sections(
      CGAL::Polygon_mesh_slicer<Mesh, Kernel>& slicer,
      const CrossSectionSampler& sampler,
      OutputIterator out)
    {
      for(int i(0); i < sampler.num_planes(); i++)
      {
        CuttingPlane cutting_plane = sampler.get_plane(i);
        find_cross_section(slicer, cutting_plane, out);
      }
    }

    template<typename CuttingPlaneIterator, typename OutputIterator>
    void find_cross_sections(
      CGAL::Polygon_mesh_slicer<Mesh, Kernel>& slicer,
      const CuttingPlaneIterator start,
      const CuttingPlaneIterator end,
      OutputIterator out)
    {
      for(CuttingPlaneIterator it = start; it != end; it++)
      {
        find_cross_section(slicer, *it, out);
      }
    }

    Point3D face_center(const Mesh& mesh, FaceIndex f)
    {
      Point3D centroid(0.0, 0.0, 0.0);
      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(f)))
      {
        centroid += Real(1.0 / 3.0) * Vector3D(Point3D(CGAL::Origin()), mesh.point(v));
      }
      return centroid;
    }

    std::size_t last_complete_depth(Octree& octree, Octree::Node_index index)
    {
      if(octree.is_leaf(index))
      {
        return 0;
      }
      else
      {
        std::size_t min = last_complete_depth(octree, octree.child(index, 0));
        for(std::size_t i(1); i < 8; i++)
        {
          min = std::min(min, last_complete_depth(octree, octree.child(index, i)));
        }
        return 1 + min;
      }
    }

    std::size_t last_complete_depth(Octree& octree)
    {
      return last_complete_depth(octree, octree.root());
    }
  } // namespace

  Point2D CuttingPlane::project(const Point3D& point) const
  {
    auto x = CGAL::scalar_product(Vector3D(origin, point), x_axis);
    auto y = CGAL::scalar_product(Vector3D(origin, point), y_axis);

    return {x, y};
  }

  Point3D RadialCrossSectionSampler::origin() const
  {
    return _origin;
  }

  int RadialCrossSectionSampler::num_planes() const
  {
    return _num_planes;
  }

  CuttingPlane RadialCrossSectionSampler::get_plane(int idx) const
  {
    idx = std::max(0, std::min(idx, _num_planes - 1));

    const double delta_angle = _num_planes > 1 ? (2.0 * M_PI) / double(_num_planes - 1) : 2.0 * M_PI;
    const double angle = idx * delta_angle;

    Vector3D plane_normal = std::cos(angle) * _u + std::sin(angle) * _v;
    Plane3D plane(_origin, plane_normal);

    Vector3D y_axis = _up;
    Vector3D x_axis = CGAL::cross_product(plane_normal, y_axis);

    return CuttingPlane{plane, _origin, x_axis, y_axis};
  }

  Point2D RadialCrossSectionSampler::project(Point3D p) const
  {
    // We need a radial projection onto the 0-th cutting plane
    CuttingPlane plane = get_plane(0);

    // Figure out height
    const auto height = CGAL::scalar_product(Vector3D(plane.origin, p), plane.y_axis);

    Point3D ref = plane.origin + height * plane.y_axis;

    // p, ref, and the projected point now all lie on a plane orthogonal to the axis
    // Target point is then the point on the 0-th cutting plane with the same distance
    // as the vector p - ref.

    const auto dist = CGAL::approximate_sqrt(Vector3D(ref, p).squared_length());

    return {dist, height};
  }

  CuttingPlane RadialCrossSectionSampler::get_plane_through_vertex(Point3D p) const
  {
    Vector3D y_axis = _up;
    Vector3D x_axis = Vector3D(_origin, p) - y_axis * CGAL::scalar_product(Vector3D(_origin, p), y_axis);
    x_axis = x_axis / CGAL::approximate_sqrt(x_axis.squared_length());
    Vector3D normal = CGAL::cross_product(x_axis, y_axis);

    return CuttingPlane{Plane3D(_origin, normal), _origin, x_axis, y_axis};
  }

  Point3D LineCrossSectionSampler::origin() const
  {
    return _start;
  }

  int LineCrossSectionSampler::num_planes() const
  {
    return _num_planes;
  }

  CuttingPlane LineCrossSectionSampler::get_plane(int idx) const
  {
    idx = std::max(0, std::min(idx, _num_planes - 1));

    const double delta = 1.0 / double(_num_planes - 1);

    const Point3D origin = _start + double(idx) * delta * Vector3D(_start, _end);

    Plane3D plane(origin, _normal);

    return CuttingPlane{plane, origin, _x_axis, _y_axis};
  }

  Point2D LineCrossSectionSampler::project(Point3D p) const
  {
    // Orthogonal projection onto 0-th cutting plane

    CuttingPlane plane = get_plane(0);

    const auto x = CGAL::scalar_product(Vector3D(plane.origin, p), plane.x_axis);
    const auto y = CGAL::scalar_product(Vector3D(plane.origin, p), plane.y_axis);

    return {x, y};
  }

  CuttingPlane LineCrossSectionSampler::get_plane_through_vertex(Point3D p) const
  {
    const auto distance = CGAL::scalar_product(_normal, Vector3D(_start, p));
    const Point3D origin = _start + distance * _normal;

    Plane3D plane(origin, _normal);

    return CuttingPlane{plane, origin, _x_axis, _y_axis};
  }

  /**
   * \brief Simplifies a polyline by merging runs of segments with similar normals
   *
   * Normals whose angles differ by at most \c threshold are considered similar.
   *
   * \param polyline The polyline to simplify
   * \param threshold Similarity threshold for normals, in degrees.
   */
  std::pair<Polygon2D, std::vector<std::size_t>> simplify_by_normal(
    const Polygon2D& polygon,
    const std::function<bool(const std::vector<Vector2D>&, const Vector2D&)>& continue_pred)
  {
    Polyline2D result;
    std::vector<std::size_t> chosen_points;

    // Find vertex of polyline with smallest inside angle as a safe starting point
    const std::size_t size = polygon.size();
    std::size_t starting_vertex = 0;
    Real minimal_angle = 360.0;
    for(std::size_t i(0); i < polygon.size(); i++)
    {
      const Point2D& a = polygon[i];
      const Point2D& b = polygon[(i + 1) % size];
      const Point2D& c = polygon[(i + 2) % size];

      Vector2D ba = a - b;
      ba = ba / CGAL::approximate_sqrt(ba.squared_length());
      Vector2D bc = c - b;
      bc = bc / CGAL::approximate_sqrt(bc.squared_length());

      const Real theta = angle(ba, bc);

      if(theta < minimal_angle)
      {
        starting_vertex = i + 1;
        minimal_angle = theta;
      }
    }

    // Build simplified result by finding sections with similar normals
    WindingOrder order = polygon.is_clockwise_oriented() ? WindingOrder::Clockwise : WindingOrder::CounterClockwise;

    std::size_t section_start = starting_vertex;
    std::size_t section_end = (starting_vertex + 1) % size;
    std::size_t candidate = (starting_vertex + 2) % size;

    Vector2D reference_normal = outside_normal(polygon[section_start], polygon[section_end], order);
    Vector2D link_normal = outside_normal(polygon[section_end], polygon[candidate], order);

    result.push_back(polygon[section_start]);
    chosen_points.push_back(section_start);

    std::vector<Vector2D> section_normals;
    section_normals.push_back(reference_normal);

    while(true)
    {
      while(continue_pred(section_normals, link_normal) && section_end != starting_vertex)
      {
        section_end = candidate;
        candidate = (candidate + 1) % size;

        const Point2D& a = polygon[section_end];
        const Point2D& b = polygon[candidate];

        section_normals.push_back(link_normal);
        link_normal = outside_normal(a, b, order);
      }

      // Latest link_normal is outside threshold
      // Merge the links from [section_start, section_end]
      // and continue next section at section_end if not all links are covered yet.

      if(section_end == starting_vertex)
      {
        // With this section we have covered all links of the original polyline.
        // The last vertex we pushed thus just gets connected to the inital vertex.

        return {Polygon2D(result.begin(), result.end()), chosen_points};
      }

      // We have not yet covered all links of the original polyline.
      // Merge the current section into a single segment and then continue.
      result.push_back(polygon[section_end]);
      chosen_points.push_back(section_end);

      section_start = section_end;
      section_end = candidate;
      candidate = (candidate + 1) % size;

      reference_normal = outside_normal(polygon[section_start], polygon[section_end], order);
      link_normal = outside_normal(polygon[section_end], polygon[candidate], order);

      section_normals.clear();
      section_normals.push_back(reference_normal);
    }

    return {Polygon2D(), std::vector<std::size_t>()};
  }

  /**
   * \brief Computes the union of cross sections of a surface mesh
   *
   * Cross sections are generated by rotating a plane.
   *
   * The axis of rotation is given by the point \c p and the vector \c up.
   * The initial cross section is given by the plane going through \c p with normal \c u.
   * The plane is then rotated 360 degrees around the axis of rotation to produce \c cross_sections cross sections.
   *
   * \returns The union of all generated cross sections.
   *
   * \param filename Path of the surface mesh file
   * \param p Point on cutting plane
   * \param up Up-direction. Cutting plane will run through axis described by p and up
   * \param u Vector orthogonal to up. Initial normal of cutting plane
   * \param cross_sections Number of cross sections
   */
  std::vector<PolygonWithHoles2D> union_of_cross_sections(const Mesh& mesh, const CrossSectionSampler& sampler)
  {
    // Slicer constructor from the mesh
    CGAL::Polygon_mesh_slicer<Mesh, Kernel> slicer(mesh);

    std::vector<PolygonWithHoles2D> all_cross_sections;

    // Running initial sampling of mesh. Collect evenly spaced cross sections and merge them
    find_cross_sections(slicer, sampler, std::back_inserter(all_cross_sections));
    std::vector<PolygonWithHoles2D> union_components = merge(all_cross_sections);

    // The inital samples should have captured the major features of the input mesh,
    // but there is no guarantee that all features have been captured.
    // As a next step the vertices of the input mesh onto the components of the union.
    // If a vertex gets projected outside of the outer boundary or inside a hole,
    // we have missed a feature of the mesh.
    // In that case we place additional cutting planes through these vertices.

    /*
    std::vector<Point> outside_vertices;

    for(const Point& vertex : mesh.points())
    {
      Point2D projected = sampler.project(vertex);
      if(is_outside_union(projected, union_components))
      {
        outside_vertices.push_back(vertex);
      }
    }

    // Merge additional feature planes into union

    // Sort vertices by distance to sampler origin.
    // For radial sampler this should save some time,
    // because more distant vertices will cover more other vertices
    // when added to the union.
    std::sort(outside_vertices.begin(), outside_vertices.end(), [&](Point& a, Point& b) {
      return Vector(sampler.origin(), a).squared_length() > Vector(sampler.origin(), b).squared_length();
    });

    int feature_planes_checked = 0;
    for(Point& vertex : outside_vertices)
    {
      Point2D projected = sampler.project(vertex);
      if(is_outside_union(projected, union_components))
      {
        CuttingPlane plane = sampler.get_plane_through_vertex(vertex);
        find_cross_section(slicer, plane, std::back_inserter(union_components));
        union_components = merge(union_components);
      }
      feature_planes_checked++;
    }
    */

    return union_components;
  }

  VolumeMesh fbm_mesh(Mesh& mesh, const FBMMeshSettings& settings)
  {
    const BoundingBox& bb = settings.bounding_box;

    // Determine local gaps and local min gap
    std::vector<Gap> gap_vec = gaps(mesh);

    // Determine orientation corrected target cell sizes
    // For each gap we can determine the target cell sizes in x, y, z direction
    // at that point

    std::vector<std::pair<FaceIndex, Point3D>> target_cell_sizes;
    target_cell_sizes.reserve(gap_vec.size());

    // Create masks to filter out gaps at self-intersecting faces
    // TODO: Move this to gap-scoring. Not done yet because calculating
    // gaps is slow and doing so would invalidate all current checkpoints.
    std::vector<std::pair<FaceIndex, FaceIndex>> self_intersections;
    CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(self_intersections));
    std::vector<bool> self_intersection_mask(mesh.number_of_faces(), false);
    for(auto& pair : self_intersections)
    {
      self_intersection_mask[static_cast<std::uint32_t>(pair.first)] = true;
      self_intersection_mask[static_cast<std::uint32_t>(pair.second)] = true;
    }

    for(const Gap& gap : gap_vec)
    {
      if(
        gap.confidence > 0.9 && !self_intersection_mask[static_cast<std::uint32_t>(gap.face)] &&
        !self_intersection_mask[static_cast<std::uint32_t>(gap.opposite_face)])
      {
        const FaceIndex start_face(gap.face);
        const FaceIndex end_face(gap.opposite_face);

        const Vector3D gap_dir = CGAL::Polygon_mesh_processing::compute_face_normal(start_face, mesh);

        const double eps = 1e-4;
        const double x_size = gap.diameter / (std::abs(gap_dir.x()) + eps);
        const double y_size = gap.diameter / (std::abs(gap_dir.y()) + eps);
        const double z_size = gap.diameter / (std::abs(gap_dir.z()) + eps);

        const double base_x_size = x_size * std::pow(2, settings.levels);
        const double base_y_size = y_size * std::pow(2, settings.levels);
        const double base_z_size = z_size * std::pow(2, settings.levels);

        Point3D target_cell_size{base_x_size, base_y_size, base_z_size};

        target_cell_sizes.emplace_back(start_face, target_cell_size);
      }
    }

    // Write target cell sizes to mesh properties
    // TODO: Only do this if properties do not exist yet
    auto x_prop = mesh.add_property_map<FaceIndex, double>("f:target_cell_size_x").first;
    auto y_prop = mesh.add_property_map<FaceIndex, double>("f:target_cell_size_y").first;
    auto z_prop = mesh.add_property_map<FaceIndex, double>("f:target_cell_size_z").first;

    for(const std::pair<FaceIndex, Point3D>& t : target_cell_sizes)
    {
      x_prop[t.first] = t.second.x();
      y_prop[t.first] = t.second.y();
      z_prop[t.first] = t.second.z();
    }

    // Build size field, i.e. pairs of (Location, TargetSize)
    std::vector<std::pair<Point3D, Point3D>> size_field;
    size_field.reserve(target_cell_sizes.size());
    for(const auto& t : target_cell_sizes)
    {
      size_field.emplace_back(face_center(mesh, t.first), t.second);
    }

    // Create octree with custom splitting predicate.
    // Split as long as any point of a cell requires a mesh-cell smaller then the current octree-cell

    Octree octree(size_field, CGAL::First_of_pair_property_map<std::pair<Point3D, Point3D>>());
    octree.refine(
      [&](auto node, const auto& tree)
      {
        const std::size_t depth = tree.depth(node);

        const double width_x = ((bb.max.x - bb.min.x) / std::pow(2, depth)) * 0.9;
        const double width_y = ((bb.max.y - bb.min.y) / std::pow(2, depth)) * 0.9;
        const double width_z = ((bb.max.z - bb.min.z) / std::pow(2, depth)) * 0.9;

        for(const auto& data : tree.data(node))
        {
          const Point3D target_size = data.second;

          if(target_size.x() < width_x)
          {
            return true;
          }
          if(target_size.y() < width_y)
          {
            return true;
          }
          if(target_size.z() < width_z)
          {
            return true;
          }
        }

        return false;
      });

    // Determine deepest complete level of octree
    const std::size_t global_refinements = std::max(last_complete_depth(octree), 2UL);

    // Create base mesh
    const auto verts_per_axis = static_cast<std::size_t>(std::pow(2, global_refinements) + 1);

    std::vector<double> slices_x(verts_per_axis);
    std::vector<double> slices_y(verts_per_axis);
    std::vector<double> slices_z(verts_per_axis);

    const double delta_x = (bb.max.x - bb.min.x) / static_cast<double>(verts_per_axis - 1);
    const double delta_y = (bb.max.y - bb.min.y) / static_cast<double>(verts_per_axis - 1);
    const double delta_z = (bb.max.z - bb.min.z) / static_cast<double>(verts_per_axis - 1);

    for(std::size_t i(0); i < verts_per_axis; i++)
    {
      slices_x[i] = bb.min.x + (static_cast<double>(i) * delta_x);
      slices_y[i] = bb.min.y + (static_cast<double>(i) * delta_y);
      slices_z[i] = bb.min.z + (static_cast<double>(i) * delta_z);
    }

    VolumeMesh
      result{slices_x.begin(), slices_x.end(), slices_y.begin(), slices_y.end(), slices_z.begin(), slices_z.end()};


    for(const std::pair<Point3D, Point3D>& p : size_field)
    {
      auto node_index = octree.locate(p.first);
      // NOTE: Because we force at least two full refinements above, we have
      // to assume at least depth 2 here
      const std::size_t depth = std::max(octree.depth(node_index), 2UL);

      if(depth <= global_refinements)
      {
      }

      const auto adaptive_refinements = static_cast<std::size_t>(
        std::ceil(static_cast<double>(depth - global_refinements) * std::log(2) / std::log(3)));

      // Assign subdivision levels to closest vertex of cell containing the size_field point
      for(std::size_t cell(0); cell < result.num_cells(); cell++)
      {
        Point min = result.vertex(result.cell(cell, 0));
        Point max = result.vertex(result.cell(cell, 7));

        if(
          min.x <= p.first.x() && p.first.x() <= max.x && min.y <= p.first.y() && p.first.y() <= max.y &&
          min.z <= p.first.z() && p.first.z() <= max.z)
        {
          std::size_t closest_idx = 0;
          Point3D closest_point = Point3D(min.x, min.y, min.z);
          double closest_distance = std::numeric_limits<double>::max();

          for(int i(0); i < 8; i++)
          {
            std::size_t vertex = result.cell(cell, i);
            Point vertex_point = result.vertex(vertex);
            Point3D v(vertex_point.x, vertex_point.y, vertex_point.z);

            double distance = (v - p.first).squared_length();

            if(distance < closest_distance)
            {
              closest_idx = vertex;
              closest_point = v;
              closest_distance = distance;
            }
          }
          std::uint64_t& current = result.subdivision_level(closest_idx);
          current = std::max(current, adaptive_refinements);
        }
      }
    }

    return result;
  }
} // namespace MeshHexer
