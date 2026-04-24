#include <meshhexer/meshhexer.hpp>

#include <cgal_types.hpp>
#include <io.hpp>
#include <macros.hpp>
#include <meshing.hpp>
#include <properties.hpp>
#include <meshhexer/types.hpp>
#include <meshhexer/config.hpp>
#include <warnings.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <numeric>
#include <sstream>
#include <unistd.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/squared_distance_2.h>

#include <Eigen/Dense>

namespace MeshHexer
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  namespace
  {
    struct WeightedSample
    {
      double value = 0.0;
      double weight = 0.0;
    };

    bool ends_with(const std::string& string, const std::string& ending)
    {
      if(ending.size() > string.size())
      {
        return false;
      }
      return std::equal(ending.rbegin(), ending.rend(), string.rbegin());
    }

    double squared_distance_xy(const Point3D& a, const Point3D& b)
    {
      const double dx = a.x() - b.x();
      const double dy = a.y() - b.y();
      return (dx * dx) + (dy * dy);
    }

    double polygon_area_xy(const std::vector<Point3D>& polyline)
    {
      if(polyline.size() < 3)
      {
        return 0.0;
      }

      double area = 0.0;
      for(std::size_t i = 0; i + 1 < polyline.size(); ++i)
      {
        area += (polyline[i].x() * polyline[i + 1].y()) - (polyline[i + 1].x() * polyline[i].y());
      }
      return std::abs(area) * 0.5;
    }

    std::optional<Point> fit_circle_xy(const std::vector<Point3D>& polyline)
    {
      std::vector<Point3D> points;
      points.reserve(polyline.size());
      for(const Point3D& p : polyline)
      {
        if(points.empty() || squared_distance_xy(points.back(), p) > 1e-24)
        {
          points.push_back(p);
        }
      }
      if(points.size() > 1 && squared_distance_xy(points.front(), points.back()) < 1e-24)
      {
        points.pop_back();
      }
      if(points.size() < 3)
      {
        return std::nullopt;
      }

      Eigen::Matrix3d ata = Eigen::Matrix3d::Zero();
      Eigen::Vector3d atb = Eigen::Vector3d::Zero();
      for(const Point3D& p : points)
      {
        const double x = p.x();
        const double y = p.y();
        const Eigen::Vector3d row(x, y, 1.0);
        ata += row * row.transpose();
        atb += -(x * x + y * y) * row;
      }

      if(std::abs(ata.determinant()) < 1e-18)
      {
        return std::nullopt;
      }

      const Eigen::Vector3d solution = ata.ldlt().solve(atb);
      return Point{-0.5 * solution[0], -0.5 * solution[1], 0.0};
    }

    double vector_length_xy(const Vector3D& v)
    {
      return std::sqrt((v.x() * v.x()) + (v.y() * v.y()));
    }

    double point_radius_xy(const Point3D& point, const Point& axis_center)
    {
      const double dx = point.x() - axis_center.x;
      const double dy = point.y() - axis_center.y;
      return std::sqrt((dx * dx) + (dy * dy));
    }

    double minimal_distance_to_axis_line(const Mesh& mesh, const Point& axis_center)
    {
      using Triangle2D = Kernel::Triangle_2;

      const Point2D axis_point(axis_center.x, axis_center.y);

      double min_squared_distance = std::numeric_limits<double>::max();
      for(FaceIndex face_index : mesh.faces())
      {
        std::array<Point2D, 3> points;
        int point_idx = 0;
        for(VertexIndex vertex_index : mesh.vertices_around_face(mesh.halfedge(face_index)))
        {
          if(point_idx < 3)
          {
            const Point3D& point = mesh.point(vertex_index);
            points[point_idx] = Point2D(point.x(), point.y());
          }
          point_idx++;
        }

        if(point_idx != 3)
        {
          continue;
        }

        const Triangle2D triangle(points[0], points[1], points[2]);
        min_squared_distance = std::min(
          min_squared_distance,
          CGAL::to_double(CGAL::squared_distance(axis_point, triangle)));
      }

      if(min_squared_distance == std::numeric_limits<double>::max())
      {
        return 0.0;
      }
      return std::sqrt(std::max(0.0, min_squared_distance));
    }

    double maximal_distance_to_axis(const Mesh& mesh, const Point& axis_center)
    {
      double max_radius = 0.0;
      for(VertexIndex vertex_index : mesh.vertices())
      {
        max_radius = std::max(max_radius, point_radius_xy(mesh.point(vertex_index), axis_center));
      }
      return max_radius;
    }
  } // namespace

  class SurfaceMesh::SurfaceMeshImpl
  {
    /// Surface mesh
    Mesh _mesh;

    /// Cached AABBTree
    std::optional<AABBTree> _aabb_tree;

    std::vector<MonitorHistogramBin> _monitor_histogram;
    std::vector<MonitorHistogramComponentStats> _monitor_histogram_component_counts;
    double _monitor_histogram_total_area = 0.0;
    double _surface_total_area = 0.0;
    bool _monitor_histogram_ready = false;
    std::size_t _monitor_hist_min_index = 0;
    std::size_t _monitor_hist_max_index = 0;

  public:
    explicit SurfaceMeshImpl(Mesh&& m) : _mesh(std::move(m))
    {
      auto validity_prop = _mesh.add_property_map<FaceIndex, std::uint32_t>("f:Validity", 0);
      if(validity_prop.second)
      {
        for(FaceIndex f : _mesh.faces())
        {
          validity_prop.first[f] = 0;
        }
      }
    }

    AABBTree& aabb_tree()
    {
      if(!_aabb_tree.has_value())
      {
        _aabb_tree = AABBTree(_mesh.faces_begin(), _mesh.faces_end(), _mesh);
      }
      return _aabb_tree.value();
    }

    /// \copydoc SurfaceMesh::bounding_box()
    BoundingBox bounding_box() const;

    /// \copydoc SurfaceMesh::num_vertices()
    std::uint32_t num_vertices() const;
    /// \copydoc SurfaceMesh::num_edges()
    std::uint32_t num_edges() const;
    /// \copydoc SurfaceMesh::num_faces()
    std::uint32_t num_faces() const;

    /// \copydoc SurfaceMesh::is_closed()
    bool is_closed() const;
    /// \copydoc SurfaceMesh::is_wound_consistently()
    bool is_wound_consistently() const;
    /// \copydoc SurfaceMesh::is_outward_oriented()
    bool is_outward_oriented() const;
    /// \copydoc SurfaceMesh::minimal_aspect_ratio()
    double minimal_aspect_ratio() const;
    /// \copydoc SurfaceMesh::maximal_aspect_ratio()
    double maximal_aspect_ratio() const;

    /// \copydoc SurfaceMesh::gaps()
    std::vector<Gap> gaps();
    /// \copydoc SurfaceMesh::min_gap()
    Gap min_gap();

    /// \copydoc SurfaceMesh::fbm_mesh()
    VolumeMesh fbm_mesh(const FBMMeshSettings& settings);

    /// \copydoc SurfaceMesh::warnings()
    MeshWarnings warnings() const;

    /// \copydoc SurfaceMesh::faces_with_short_normals()
    std::vector<std::size_t> faces_with_short_normals(double min_length);

    /// \copydoc SurfaceMesh::round_geometry_analysis()
    Result<RoundGeometryAnalysisResult, std::string> round_geometry_analysis(
      const std::vector<RoundGeometryInflow>& inflows,
      const RoundGeometryAnalysisConfig& config) const;

    /// \copydoc SurfaceMesh::translate()
    void translate(double dx, double dy, double dz);

    /// \copydoc SurfaceMesh::write_monitor_histogram()
    void write_monitor_histogram(std::ostream& stream) const;

    /// \copydoc SurfaceMesh::write_to_file()
    Result<void, std::string> write_to_file(const std::string& filename);

  private:
    /**
     * \brief Ensures the given property exists on the mesh
     *
     * Does nothing if the property already exists.
     * Otherwise it calls the given function.
     */
    template<typename Index, typename T, typename Fn>
    void ensure_property(const std::string& property_name, Fn&& fn)
    {
      auto prop = _mesh.add_property_map<Index, T>(property_name, T{});

      if(prop.second)
      {
        std::forward<Fn>(fn)();
      }
    }

    void prepare_for_min_gap();
    void write_monitor_histogram_file(const std::filesystem::path& reference_path) const;
  };

  void SurfaceMesh::SurfaceMeshImpl::prepare_for_min_gap()
  {
    // Ensure vertex normals are available
    ensure_property<VertexIndex, Vector3D>("v:normals", [&]() { compute_vertex_normals(_mesh); });

    // Ensure maximal dihedral angles are available
    ensure_property<FaceIndex, double>("f:dihedral_angle", [&]() { compute_max_dihedral_angle(_mesh); });

    // Ensure maximal inscribed spheres are available
    ensure_property<FaceIndex, double>("f:MIS_diameter", [&]() { maximal_inscribed_spheres(_mesh, aabb_tree()); });
    ensure_property<FaceIndex, double>("f:normaldistance", [&]() { normal_distances(_mesh, aabb_tree()); });
    ensure_property<FaceIndex, std::uint32_t>(
      "f:normaldistance_target",
      [&]() { normal_distances(_mesh, aabb_tree()); });
    // Calculate maximum search distances for topological distances

    Mesh::Property_map<FaceIndex, double> max_distances =
      _mesh.add_property_map<FaceIndex, double>("f:max_search_distance", 0.0).first;

    Mesh::Property_map<FaceIndex, double> diameters =
      _mesh.add_property_map<FaceIndex, double>("f:MIS_diameter", 0.0).first;

    const double ms = mesh_size(_mesh);
    const MinGapConfig& config = min_gap_config();
    const double relative_minimal_search_radius = config.relative_min_search_radius;
    for(FaceIndex f : _mesh.faces())
    {
      max_distances[f] = std::max(M_PI * diameters[f], relative_minimal_search_radius * ms);
    }

    // Ensure topological distances are available
    ensure_property<FaceIndex, double>(
      "f:topological_distance",
      [&]() { topological_distances(_mesh, "f:MIS_id", "f:max_search_distance"); });

    update_validity_from_neighbor_diameters(_mesh);
    update_validity_from_neighbor_normals(_mesh);
    update_validity_from_small_angles(_mesh, config.small_angle_limit_deg);

    // Ensure gap scores are available
    score_gaps(_mesh);
  }

  std::vector<Gap> SurfaceMesh::SurfaceMeshImpl::gaps()
  {
    prepare_for_min_gap();
    return MeshHexer::gaps(_mesh);
  }

  Gap SurfaceMesh::SurfaceMeshImpl::min_gap()
  {
    prepare_for_min_gap();
    Gap gap = MeshHexer::min_gap(_mesh);

    monitor_distances(_mesh);
    _monitor_histogram = monitor_histogram(_mesh, gap.diameter);
    _monitor_histogram_component_counts = monitor_histogram_component_counts(
      _mesh,
      _monitor_histogram,
      gap.diameter,
      min_gap_config().histogram_min_fraction_connected_threshold);
    _monitor_hist_min_index = 0;
    _monitor_hist_max_index = 0;
    if(!_monitor_histogram.empty())
    {
      auto [min_idx, max_idx] = monitor_histogram_min_max_indices(_monitor_histogram);
      _monitor_hist_min_index = std::min(min_idx, _monitor_histogram.size() - 1);
      _monitor_hist_max_index = std::min(max_idx, _monitor_histogram.size() - 1);
    }
    _monitor_histogram_total_area = 0.0;
    for(const MonitorHistogramBin& bin : _monitor_histogram)
    {
      _monitor_histogram_total_area += bin.area;
    }
    _surface_total_area = 0.0;
    for(FaceIndex f : _mesh.faces())
    {
      _surface_total_area += PMP::face_area(f, _mesh);
    }
    _monitor_histogram_ready = true;

    double adjusted_min = adjusted_min_gap_from_histogram(_monitor_histogram, gap.diameter);
    gap.diameter = adjusted_min;

    invalidate_monitors_below(_mesh, gap.diameter);

    gap.bin_span = static_cast<std::uint32_t>((_monitor_hist_max_index >= _monitor_hist_min_index)
                                               ? (_monitor_hist_max_index - _monitor_hist_min_index)
                                               : 0);
    const MinGapConfig& config_span = min_gap_config();
    const std::uint32_t min_span = static_cast<std::uint32_t>(config_span.min_histogram_span_bins);
    const std::uint32_t max_span = static_cast<std::uint32_t>(config_span.max_histogram_span_bins);
    gap.bin_span = std::max(min_span, std::min(max_span, gap.bin_span));

    return gap;
  }

  VolumeMesh SurfaceMesh::SurfaceMeshImpl::fbm_mesh(const FBMMeshSettings& settings)
  {
    prepare_for_min_gap();
    return MeshHexer::fbm_mesh(_mesh, settings);
  }

  Result<void, std::string> SurfaceMesh::SurfaceMeshImpl::write_to_file(const std::string& filename)
  {
    using ResultType = Result<void, std::string>;

    if(!ends_with(filename, ".off") && !ends_with(filename, ".ply") && !ends_with(filename, ".vtu"))
    {
      return ResultType::err("Can only write .off, .ply or .vtu files");
    }

    if(ends_with(filename, ".off"))
    {
      if(!CGAL::IO::write_polygon_mesh(filename, _mesh))
      {
        return ResultType::err("Failed to write mesh to opened file");
      }
      return {};
    }

    std::ofstream output(filename);

    if(!output)
    {
      return ResultType::err("Failed to open file for writing");
    }

    if(ends_with(filename, ".ply") && !CGAL::IO::write_PLY(output, _mesh))
    {
      return ResultType::err("Failed to write mesh to opened file");
    }

    if(ends_with(filename, ".vtu"))
    {
      write_vtu(output, _mesh);

      if(output.bad())
      {
        return ResultType::err("Failed to write mesh to opened file");
      }

      write_monitor_histogram_file(filename);
    }

    return {};
  }

  void SurfaceMesh::SurfaceMeshImpl::translate(double dx, double dy, double dz)
  {
    for(VertexIndex vertex_index : _mesh.vertices())
    {
      const Point3D& point = _mesh.point(vertex_index);
      _mesh.point(vertex_index) = Point3D(point.x() + dx, point.y() + dy, point.z() + dz);
    }
  }

  void SurfaceMesh::SurfaceMeshImpl::write_monitor_histogram_file(const std::filesystem::path& reference_path) const
  {
    if(!_monitor_histogram_ready)
    {
      return;
    }

    std::filesystem::path directory = reference_path.parent_path();
    if(directory.empty())
    {
      directory = ".";
    }
    std::filesystem::path histogram_path = directory / "size_distribution_histogram.txt";

    std::ofstream out(histogram_path);
    if(!out)
    {
      return;
    }

    write_monitor_histogram(out);
  }

  void SurfaceMesh::SurfaceMeshImpl::write_monitor_histogram(std::ostream& out) const
  {
    if(!_monitor_histogram_ready)
    {
      return;
    }

    const int width_label = 8;
    const int width_start = 18;
    const int width_end = 18;
    const int width_area = 18;
    const int width_percent = 14;
    const int width_components = 14;
    const int width_connected = 16;
    const int width_isolated = 16;
    const int width_max_component = 18;

    const int total_width = width_label + width_start + width_end + width_area + width_percent +
      width_components + width_connected + width_isolated + width_max_component;

    std::size_t marked_min = _monitor_hist_min_index;
    std::size_t marked_max = _monitor_hist_max_index;
    if(!_monitor_histogram.empty())
    {
      marked_min = std::min(marked_min, _monitor_histogram.size() - 1);
      marked_max = std::min(marked_max, _monitor_histogram.size() - 1);
    }
    else
    {
      marked_min = marked_max = std::numeric_limits<std::size_t>::max();
    }

    out << std::left << std::setw(width_label) << "Mark"
        << std::right << std::setw(width_start) << "Bin Start"
        << std::right << std::setw(width_end) << "Bin End"
        << std::right << std::setw(width_area) << "Area"
        << std::right << std::setw(width_percent) << "Percentage"
        << std::right << std::setw(width_components) << "Collections"
        << std::right << std::setw(width_connected) << "Connected %"
        << std::right << std::setw(width_isolated) << "Isolated %"
        << std::right << std::setw(width_max_component) << "Max Component %"
        << "\n";

    out << std::string(total_width, '-') << "\n";

    const auto format_value = [](double value, int precision) {
      std::ostringstream ss;
      ss << std::fixed << std::setprecision(precision) << value;
      return ss.str();
    };

    for(std::size_t i = 0; i < _monitor_histogram.size(); ++i)
    {
      const MonitorHistogramBin& bin = _monitor_histogram[i];
      double percentage = 0.0;
      if(_monitor_histogram_total_area > 0.0)
      {
        percentage = (bin.area / _monitor_histogram_total_area) * 100.0;
      }
      double connected_percentage = 0.0;
      double isolated_percentage = 0.0;
      double max_component_percentage = 0.0;
      if(i < _monitor_histogram_component_counts.size() && _monitor_histogram_total_area > 0.0)
      {
        connected_percentage = (_monitor_histogram_component_counts[i].connected_area / _monitor_histogram_total_area) * 100.0;
        isolated_percentage = (_monitor_histogram_component_counts[i].isolated_area / _monitor_histogram_total_area) * 100.0;
        max_component_percentage =
          (_monitor_histogram_component_counts[i].max_component_area / _monitor_histogram_total_area) * 100.0;
      }

      std::string label;
      if(i == marked_min && i == marked_max)
      {
        label = "MIN/MAX";
      }
      else if(i == marked_min)
      {
        label = "MIN";
      }
      else if(i == marked_max)
      {
        label = "MAX";
      }

      out << std::left << std::setw(width_label) << label
          << std::right << std::setw(width_start) << format_value(bin.lower, 6)
          << std::right << std::setw(width_end) << format_value(bin.upper, 6)
          << std::right << std::setw(width_area) << format_value(bin.area, 6)
          << std::right << std::setw(width_percent - 1) << format_value(percentage, 2) << "%"
          << std::right << std::setw(width_components)
          << ((i < _monitor_histogram_component_counts.size()) ? std::to_string(_monitor_histogram_component_counts[i].collections) : "0")
          << std::right << std::setw(width_connected - 1) << format_value(connected_percentage, 2) << "%"
          << std::right << std::setw(width_isolated - 1) << format_value(isolated_percentage, 2) << "%"
          << std::right << std::setw(width_max_component - 1) << format_value(max_component_percentage, 2) << "%"
          << "\n";
    }
  }

  BoundingBox SurfaceMesh::SurfaceMeshImpl::bounding_box() const
  {
    return MeshHexer::bounding_box(_mesh);
  }

  std::uint32_t SurfaceMesh::SurfaceMeshImpl::num_vertices() const
  {
    return _mesh.num_vertices();
  }

  std::uint32_t SurfaceMesh::SurfaceMeshImpl::num_edges() const
  {
    return _mesh.num_edges();
  }

  std::uint32_t SurfaceMesh::SurfaceMeshImpl::num_faces() const
  {
    return _mesh.num_faces();
  }

  bool SurfaceMesh::SurfaceMeshImpl::is_closed() const
  {
    return CGAL::is_closed(_mesh);
  }

  bool SurfaceMesh::SurfaceMeshImpl::is_wound_consistently() const
  {
    return MeshHexer::is_wound_consistently(_mesh);
  }

  bool SurfaceMesh::SurfaceMeshImpl::is_outward_oriented() const
  {
    XASSERT(is_closed());
    XASSERT(is_wound_consistently());
    return CGAL::Polygon_mesh_processing::is_outward_oriented(_mesh);
  }

  double SurfaceMesh::SurfaceMeshImpl::minimal_aspect_ratio() const
  {
    double min_aspect_ratio = std::numeric_limits<double>::max();
    for(MeshHexer::FaceIndex f : _mesh.faces())
    {
      double ratio = CGAL::Polygon_mesh_processing::face_aspect_ratio(f, _mesh);
      min_aspect_ratio = std::min(min_aspect_ratio, ratio);
    }

    return min_aspect_ratio;
  }

  double SurfaceMesh::SurfaceMeshImpl::maximal_aspect_ratio() const
  {
    double max_aspect_ratio = 0;
    for(MeshHexer::FaceIndex f : _mesh.faces())
    {
      double ratio = CGAL::Polygon_mesh_processing::face_aspect_ratio(f, _mesh);
      max_aspect_ratio = std::max(max_aspect_ratio, ratio);
    }

    return max_aspect_ratio;
  }

  MeshWarnings SurfaceMesh::SurfaceMeshImpl::warnings() const
  {
    MeshWarnings ws;
    create_warnings(_mesh, ws);
    return ws;
  }

  std::vector<std::size_t> SurfaceMesh::SurfaceMeshImpl::faces_with_short_normals(double min_length)
  {
    std::vector<FaceIndex> invalid_faces = faces_with_short_vertex_normals(_mesh, min_length);
    std::vector<std::size_t> indices;
    indices.reserve(invalid_faces.size());

    for(FaceIndex f : invalid_faces)
    {
      indices.push_back(static_cast<std::size_t>(f));
    }

    return indices;
  }

  Result<RoundGeometryAnalysisResult, std::string> SurfaceMesh::SurfaceMeshImpl::round_geometry_analysis(
    const std::vector<RoundGeometryInflow>& inflows,
    const RoundGeometryAnalysisConfig& config) const
  {
    using ResultType = Result<RoundGeometryAnalysisResult, std::string>;

    if(_mesh.is_empty())
    {
      return ResultType::err("Round geometry analysis requires a non-empty mesh.");
    }
    if(inflows.empty())
    {
      return ResultType::err("Round geometry analysis requires at least one inflow description.");
    }

    const BoundingBox bb = MeshHexer::bounding_box(_mesh);
    const double z_extent = bb.max.z - bb.min.z;
    if(z_extent <= 0.0)
    {
      return ResultType::err("Round geometry analysis requires a positive z extent.");
    }

    const double z_slice = bb.max.z - (config.top_slice_relative_epsilon * z_extent);
    const double close_tol_sq = std::max(1e-24, std::pow(mesh_size(_mesh) * 1e-8, 2.0));

    CGAL::Polygon_mesh_slicer<Mesh, Kernel> slicer(_mesh);
    Polylines3D polylines;
    slicer(Plane3D(Point3D(0.0, 0.0, z_slice), Vector3D(0.0, 0.0, 1.0)), std::back_inserter(polylines));

    std::vector<Point3D> selected_loop;
    double selected_area = -1.0;
    for(const Polyline3D& polyline : polylines)
    {
      if(polyline.size() < 3)
      {
        continue;
      }

      std::vector<Point3D> loop(polyline.begin(), polyline.end());
      if(squared_distance_xy(loop.front(), loop.back()) > close_tol_sq)
      {
        continue;
      }
      if(squared_distance_xy(loop.front(), loop.back()) > 0.0)
      {
        loop.push_back(loop.front());
      }

      const double area = polygon_area_xy(loop);
      if(area > selected_area)
      {
        selected_area = area;
        selected_loop = std::move(loop);
      }
    }

    if(selected_loop.empty())
    {
      return ResultType::err("Failed to extract a closed outer loop from the z-max slice.");
    }

    const std::optional<Point> maybe_center = fit_circle_xy(selected_loop);
    if(!maybe_center.has_value())
    {
      return ResultType::err("Failed to fit a circle to the z-max slice.");
    }
    const Point axis_center = maybe_center.value();

    double outer_radius = maximal_distance_to_axis(_mesh, axis_center);
    double inner_radius = minimal_distance_to_axis_line(_mesh, axis_center);

    std::vector<double> outer_radius_limits;
    std::vector<double> extrusion_lengths;
    outer_radius_limits.reserve(inflows.size());
    extrusion_lengths.reserve(inflows.size());

    const double axial_dot_threshold = std::cos((config.axial_inflow_tolerance_deg / 180.0) * M_PI);
    const double radial_constraint_tol = std::max(mesh_size(_mesh) * 1e-6, 1e-9);
    bool z_min_limited = false;

    for(const RoundGeometryInflow& inflow : inflows)
    {
      if(!inflow.has_center || !inflow.has_normal)
      {
        continue;
      }

      const Point3D inflow_center(inflow.center.x, inflow.center.y, inflow.center.z);
      Vector3D inflow_normal(inflow.normal.x, inflow.normal.y, inflow.normal.z);
      const double inflow_normal_length = std::sqrt(inflow_normal.squared_length());
      if(inflow_normal_length <= 0.0)
      {
        continue;
      }
      inflow_normal = inflow_normal / inflow_normal_length;

      const double inflow_radius = point_radius_xy(inflow_center, axis_center);
      const double inflow_normal_xy_length = vector_length_xy(inflow_normal);
      if(inflow_radius > radial_constraint_tol && inflow_normal_xy_length > radial_constraint_tol)
      {
        const double radial_x = (inflow_center.x() - axis_center.x) / inflow_radius;
        const double radial_y = (inflow_center.y() - axis_center.y) / inflow_radius;
        const double radial_alignment =
          ((inflow_normal.x() * radial_x) + (inflow_normal.y() * radial_y)) / inflow_normal_xy_length;

        if(radial_alignment < 0.0)
        {
          outer_radius_limits.push_back(inflow_radius);
        }
      }

      double min_projection = std::numeric_limits<double>::max();
      for(VertexIndex vertex_index : _mesh.vertices())
      {
        const Point3D point = _mesh.point(vertex_index);
        const double projection = CGAL::scalar_product(Vector3D(inflow_center, point), inflow_normal);
        min_projection = std::min(min_projection, projection);
      }
      if(min_projection != std::numeric_limits<double>::max())
      {
        extrusion_lengths.push_back(std::max(0.0, -min_projection));
      }

      if(inflow_normal.z() > axial_dot_threshold)
      {
        z_min_limited = true;
      }
    }

    if(!outer_radius_limits.empty())
    {
      outer_radius = std::min(outer_radius, *std::min_element(outer_radius_limits.begin(), outer_radius_limits.end()));
    }
    if(outer_radius <= 0.0)
    {
      return ResultType::err("Round geometry analysis produced a non-positive outer radius.");
    }
    if(extrusion_lengths.empty())
    {
      return ResultType::err("Failed to estimate extrusion length from inflows.");
    }

    std::vector<double> sorted_lengths = extrusion_lengths;
    std::sort(sorted_lengths.begin(), sorted_lengths.end());
    const double extrusion_length = sorted_lengths[sorted_lengths.size() / 2];
    const double extrusion_tol = std::max(
      mesh_size(_mesh) * 1e-8,
      config.extrusion_length_consistency_relative_tolerance * std::max(extrusion_length, 1.0));

    RoundGeometryAnalysisResult result;
    result.axis_center = axis_center;
    result.axis_aligned_to_origin =
      (std::abs(axis_center.x) <= config.axis_alignment_tolerance) &&
      (std::abs(axis_center.y) <= config.axis_alignment_tolerance);
    result.top_slice_z = z_slice;
    result.outer_diameter = 2.0 * outer_radius;
    result.inner_diameter = std::max(0.0, 2.0 * inner_radius);
    result.inner_to_outer_ratio =
      (result.outer_diameter > 0.0) ? (result.inner_diameter / result.outer_diameter) : 0.0;
    result.classification =
      (result.inner_to_outer_ratio < config.full_cylinder_inner_to_outer_threshold) ?
        RoundGeometryClassification::FullCylinder :
        RoundGeometryClassification::HollowCylinder;
    result.extrusion_length = extrusion_length;
    result.extrusion_length_consistent = std::all_of(
      extrusion_lengths.begin(),
      extrusion_lengths.end(),
      [&](double value) { return std::abs(value - extrusion_length) <= extrusion_tol; });
    result.extrusion_length_samples = std::move(extrusion_lengths);
    result.z_min_limited_by_axial_inflow = z_min_limited;
    result.z_min_physical = bb.min.z + (z_min_limited ? extrusion_length : 0.0);
    result.z_max_physical = bb.max.z - extrusion_length;

    return ResultType::ok(std::move(result));
  }

  SurfaceMesh::SurfaceMesh() = default;

  SurfaceMesh::SurfaceMesh(std::unique_ptr<SurfaceMesh::SurfaceMeshImpl> ptr) : impl(std::move(ptr))
  {
  }

  SurfaceMesh::SurfaceMesh(SurfaceMesh&&) noexcept = default;
  SurfaceMesh& SurfaceMesh::operator=(SurfaceMesh&&) noexcept = default;
  SurfaceMesh::~SurfaceMesh() = default;

  std::vector<Gap> SurfaceMesh::gaps()
  {
    return impl->gaps();
  }

  Gap SurfaceMesh::min_gap()
  {
    return impl->min_gap();
  }

  VolumeMesh SurfaceMesh::fbm_mesh(const FBMMeshSettings& settings)
  {
    return impl->fbm_mesh(settings);
  }

  Result<void, std::string> SurfaceMesh::write_to_file(const std::string& filename)
  {
    return impl->write_to_file(filename);
  }

  BoundingBox SurfaceMesh::bounding_box() const
  {
    return impl->bounding_box();
  }

  std::uint32_t SurfaceMesh::num_vertices() const
  {
    return impl->num_vertices();
  }

  std::uint32_t SurfaceMesh::num_edges() const
  {
    return impl->num_edges();
  }

  std::uint32_t SurfaceMesh::num_faces() const
  {
    return impl->num_faces();
  }

  bool SurfaceMesh::is_closed() const
  {
    return impl->is_closed();
  }

  bool SurfaceMesh::is_wound_consistently() const
  {
    return impl->is_wound_consistently();
  }

  bool SurfaceMesh::is_outward_oriented() const
  {
    return impl->is_outward_oriented();
  }

  double SurfaceMesh::minimal_aspect_ratio() const
  {
    return impl->minimal_aspect_ratio();
  }

  double SurfaceMesh::maximal_aspect_ratio() const
  {
    return impl->maximal_aspect_ratio();
  }

  MeshWarnings SurfaceMesh::warnings() const
  {
    return impl->warnings();
  }

  std::vector<std::size_t> SurfaceMesh::faces_with_short_normals(double min_length)
  {
    return impl->faces_with_short_normals(min_length);
  }

  Result<RoundGeometryAnalysisResult, std::string> SurfaceMesh::round_geometry_analysis(
    const std::vector<RoundGeometryInflow>& inflows,
    const RoundGeometryAnalysisConfig& config)
  {
    return impl->round_geometry_analysis(inflows, config);
  }

  void SurfaceMesh::translate(double dx, double dy, double dz)
  {
    impl->translate(dx, dy, dz);
  }

  void SurfaceMesh::write_monitor_histogram(std::ostream& stream) const
  {
    impl->write_monitor_histogram(stream);
  }

  Result<SurfaceMesh, std::string> load_from_file(const std::string& filename, bool triangulate)
  {
    using ResultType = Result<SurfaceMesh, std::string>;

    MeshHexer::Mesh mesh;
    if(ends_with(filename, ".ply"))
    {
      std::ifstream mesh_file(filename);
      std::string comment;
      if(!CGAL::IO::read_PLY(mesh_file, mesh, comment, true))
      {
        return ResultType::err("Failed to read mesh " + filename);
      }
    }
    else if(ends_with(filename, ".vtu"))
    {
      std::ifstream mesh_file(filename);
      Result<Mesh, std::string> read_result = read_vtu(mesh_file);
      if(read_result.is_err())
      {
        return Result<SurfaceMesh, std::string>::err(std::move(read_result).take_err());
      }
      mesh = std::move(read_result).take_ok();
    }
    else if(!PMP::IO::read_polygon_mesh(filename, mesh))
    {
      return ResultType::err("Failed to read mesh " + filename);
    }

    if(CGAL::is_empty(mesh))
    {
      return ResultType::err("Mesh " + filename + " is empty.");
    }

    bool is_triangle_mesh = CGAL::is_triangle_mesh(mesh);
    if(!is_triangle_mesh && triangulate)
    {
      PMP::triangulate_faces(mesh);
    }
    else if(!is_triangle_mesh && !triangulate)
    {
      return ResultType::err("Mesh " + filename + " is not a triangle mesh.");
    }

    auto impl = std::make_unique<SurfaceMesh::SurfaceMeshImpl>(std::move(mesh));
    SurfaceMesh smesh(std::move(impl));

    return ResultType::ok(std::move(smesh));
  }
} // namespace MeshHexer
