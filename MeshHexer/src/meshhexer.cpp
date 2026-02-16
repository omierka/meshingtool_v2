#include <meshhexer/meshhexer.hpp>

#include <cgal_types.hpp>
#include <io.hpp>
#include <macros.hpp>
#include <meshing.hpp>
#include <properties.hpp>
#include <meshhexer/types.hpp>
#include <warnings.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <optional>
#include <unistd.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace MeshHexer
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  namespace
  {
    bool ends_with(const std::string& string, const std::string& ending)
    {
      if(ending.size() > string.size())
      {
        return false;
      }
      return std::equal(ending.rbegin(), ending.rend(), string.rbegin());
    }
  } // namespace

  class SurfaceMesh::SurfaceMeshImpl
  {
    /// Surface mesh
    Mesh _mesh;

    /// Cached AABBTree
    std::optional<AABBTree> _aabb_tree;

    std::vector<MonitorHistogramBin> _monitor_histogram;
    double _monitor_histogram_total_area = 0.0;
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
    // Search for at least 0.5% of the mesh size
    constexpr double relative_minimal_search_radius = 0.005;
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
    update_validity_from_small_angles(_mesh);

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
    _monitor_histogram_ready = true;

    double adjusted_min = adjusted_min_gap_from_histogram(_monitor_histogram, gap.diameter);
    gap.diameter = adjusted_min;

    invalidate_monitors_below(_mesh, gap.diameter);

    gap.bin_span = static_cast<std::uint32_t>((_monitor_hist_max_index >= _monitor_hist_min_index)
                                               ? (_monitor_hist_max_index - _monitor_hist_min_index)
                                               : 0);
    gap.bin_span = std::max(1u, std::min(3u, gap.bin_span));

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

    if(!ends_with(filename, ".ply") && !ends_with(filename, ".vtu"))
    {
      return ResultType::err("Can only write .ply or .vtu files");
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

    const int width_label = 8;
    const int width_start = 18;
    const int width_end = 18;
    const int width_area = 18;
    const int width_percent = 14;

    const int total_width = width_label + width_start + width_end + width_area + width_percent;

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
          << std::right << std::setw(width_percent - 1) << format_value(percentage, 2) << "%\n";
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
