#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <optional>
#include <variant>
#include <vector>

namespace MeshHexer
{
  /**
   * \brief 3D Point class
   */
  struct Point
  {
    /// x component
    double x;
    /// y component
    double y;
    /// z component
    double z;
  };

  /**
   * \brief Axis aligned bounding box
   *
   * Consist of two points, describing the lexicographically smallest and largest points of the bounding box.
   */
  struct BoundingBox
  {
    /// Lower point
    Point min;
    /// Upper point
    Point max;
  };

  /**
   * \brief Possible interior diameter of a surface mesh
   */
  struct Gap
  {
    /// Starting face of gap
    std::size_t face;

    /// Face on opposite side of gap
    std::size_t opposite_face;

    /// Diameter of this gap
    double diameter;

    /// Score of this gap
    double confidence;

    /// Constructor
    Gap(std::size_t f, std::size_t of, double d, double c) :
      face(f),
      opposite_face(of),
      diameter(d),
      confidence(c)
    {
    }
  };

  /**
  * \brief Description of one slice of a structured mesh
  */
  struct Slice
  {
    /// Starting coordinate
    double coord = 0.0;

    /// Subdivision level for this slice
    std::uint64_t subdivision_level = 0;

    /// Constructor
    Slice() = default;

    /// Constructor
    Slice(double c, std::uint64_t lvl) : coord(c), subdivision_level(lvl)
    {
    }
  };

  /**
   * \brief Configuration options for creating meshes for use with the finite-boundary-method.
   */
  struct FBMMeshSettings
  {
    /**
     * \brief Bounding box for the fbm mesh
     */
    BoundingBox bounding_box;

    /**
     * \brief Number of levels in the intended mesh hierarchy
     *
     * The mesh generator uses this number to work out the adaptive refinement required,
     * such that the finest mesh of the hierarchy hits all local min gaps.
     */
    std::uint64_t levels;
  };

  /**
   * \brief A structured volume mesh
   *
   * Intended to represent non-fitted meshes for the fictitious-boundary-method.
   */
  class VolumeMesh
  {
  public:

    /// Slice iterator type
    using Iterator = std::vector<double>::iterator;

    /// Const slice iterator type
    using ConstIterator = std::vector<double>::const_iterator;

  private:
    // Grid points, including start and end points.
    std::vector<double> _xs;
    std::vector<double> _ys;
    std::vector<double> _zs;

    std::vector<std::uint64_t> _subdivision_levels;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] nx Number of elements along x-axis
     * \param[in] ny Number of elements along y-axis
     * \param[in] nz Number of elements along z-axis
     */
    VolumeMesh(std::size_t nx, std::size_t ny, std::size_t nz) :
      _xs(nx),
      _ys(ny),
      _zs(nz),
      _subdivision_levels(nx * ny * nz, 0)
    {
    }

    /**
     * \brief Constructor
     *
     * \tparam SliceIter Iterator over slices
     *
     * \param[in] x_begin Begin of x-slices
     * \param[in] x_end   End of x-slices
     * \param[in] y_begin Begin of y-slices
     * \param[in] y_end   End of y-slices
     * \param[in] z_begin Begin of z-slices
     * \param[in] z_end   End of z-slices
     */
    template<typename DoubleIter>
    VolumeMesh(
      DoubleIter x_begin,
      DoubleIter x_end,
      DoubleIter y_begin,
      DoubleIter y_end,
      DoubleIter z_begin,
      DoubleIter z_end) :
      _xs(x_begin, x_end),
      _ys(y_begin, y_end),
      _zs(z_begin, z_end),
      _subdivision_levels(_xs.size() * _ys.size() * _zs.size(), 0)
    {
    }

    /// Range-begin for x-slices
    Iterator xs_begin()
    {
      return _xs.begin();
    }

    /// Range-end for x-slices
    Iterator xs_end()
    {
      return _xs.end();
    }

    /// Range-begin for y-slices
    Iterator ys_begin()
    {
      return _ys.begin();
    }

    /// Range-end for y-slices
    Iterator ys_end()
    {
      return _ys.end();
    }

    /// Range-begin for z-slices
    Iterator zs_begin()
    {
      return _zs.begin();
    }

    /// Range-end for z-slices
    Iterator zs_end()
    {
      return _zs.end();
    }

    /// Returns number of vertices of mesh
    std::size_t num_vertices() const
    {
      return _xs.size() * _ys.size() * _zs.size();
    }

    /// Returns number of edges of mesh
    std::size_t num_edges() const
    {
      return ((_xs.size() - 1) * _ys.size() * _zs.size()) + ((_ys.size() - 1) * _xs.size() * _zs.size()) +
             ((_zs.size() - 1) * _xs.size() * _ys.size());
    }

    /// Returns number of faces of mesh
    std::size_t num_faces() const
    {
      return (_xs.size() * (_ys.size() - 1) * (_zs.size() - 1)) + (_ys.size() * (_xs.size() - 1) * (_zs.size() - 1)) +
             (_zs.size() * (_xs.size() - 1) * (_ys.size() - 1));
    }

    /// Returns number of cells of mesh
    std::size_t num_cells() const
    {
      return (_xs.size() - 1) * (_ys.size() - 1) * (_zs.size() - 1);
    }

    /**
     * \brief Get vertex coordinates
     *
     * \param[in] idx Index of vertex
     *
     * \returns The vertex coordinates of the vertex
     */
    MeshHexer::Point vertex(std::size_t idx) const
    {
      std::size_t idx_z = idx / (_xs.size() * _ys.size());
      std::size_t idx_y = idx % (_xs.size() * _ys.size()) / _xs.size();
      std::size_t idx_x = idx % (_xs.size() * _ys.size()) % _xs.size();

      return MeshHexer::Point{_xs[idx_x], _ys[idx_y], _zs[idx_z]};
    }

    /**
     * \brief Vertex at edge adjacencies
     *
     * \param[in] i Edge index
     * \param[in] j Vertex index (0 or 1)
     *
     * \returns The index of the j-th vertex of the i-th edge
     */
    std::size_t edge(std::size_t i, std::size_t j) const
    {
      // NOTE(mmuegge): Copied from FEAT3 (kernel/geometry/intern/struct_index_mapping.hpp)
      // FEATs structured mesh counts fences, rather than fence posts.
      // Hence the +- 1 everywhere

      const std::array<std::size_t, 3> num_slices = {_xs.size() - 1, _ys.size() - 1, _zs.size() - 1};

      // auxiliary index
      std::size_t pos = 0;

      // z-direction
      if(
        i >= num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1) +
               (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0])
      {
        pos = i - num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1) -
              (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0];
        return (pos / num_slices[2]) + ((pos % num_slices[2] + j) * (num_slices[0] + 1) * (num_slices[1] + 1));
      }
      // y-direction
      else if(i >= (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0])
      {
        pos = i - (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0];
        std::size_t x = (pos / num_slices[1]) % (num_slices[0] + 1);
        std::size_t y = pos % num_slices[1];
        std::size_t z = (pos / num_slices[1]) / (num_slices[0] + 1);
        return (z * (num_slices[0] + 1) * (num_slices[1] + 1)) + ((y + j) * (num_slices[0] + 1)) + x;
      }
      // x-direction
      else
      {
        return i + j + (i / num_slices[0]);
      }
    }

    /**
     * \brief Vertex at face adjacencies
     *
     * \param[in] i Face index
     * \param[in] j Vertex index (0, 1, 2, 3)
     *
     * \returns The index of the j-th vertex of the i-th face
     */
    std::size_t face(std::size_t i, std::size_t j) const
    {
      // NOTE(mmuegge): Copied from FEAT3 (kernel/geometry/intern/struct_index_mapping.hpp)
      // FEATs structured mesh counts fences, rather than fence posts.
      // Hence the +- 1 everywhere

      const std::array<std::size_t, 3> num_slices = {_xs.size() - 1, _ys.size() - 1, _zs.size() - 1};

      std::size_t pos = 0;
      if(i < num_slices[0] * num_slices[1] * (num_slices[2] + 1))
      {
        return i + (i / num_slices[0]) + ((i / (num_slices[0] * num_slices[1])) * (num_slices[0] + 1)) + (j % 2) +
               ((j / 2) * (num_slices[0] + 1));
      }
      else if(
        i >= num_slices[0] * num_slices[1] * (num_slices[2] + 1) &&
        i < num_slices[0] * (num_slices[1] * (num_slices[2] + 1) + num_slices[2] * (num_slices[1] + 1)))
      {
        pos = i - num_slices[0] * num_slices[1] * (num_slices[2] + 1);
        return (pos % num_slices[0]) + ((pos / (num_slices[0] * num_slices[2])) * (num_slices[0] + 1)) +
               ((pos % (num_slices[0] * num_slices[2])) / num_slices[0] * (num_slices[0] + 1) * (num_slices[1] + 1)) +
               (j % 2) + ((j / 2) * (num_slices[0] + 1) * (num_slices[1] + 1));
      }
      else
      {
        pos = i - num_slices[0] * (num_slices[1] * (num_slices[2] + 1) + num_slices[2] * (num_slices[1] + 1));
        return (pos / (num_slices[1] * num_slices[2])) + ((pos % num_slices[1]) * (num_slices[0] + 1)) +
               (((pos / num_slices[1]) % num_slices[2]) * (num_slices[0] + 1) * (num_slices[1] + 1)) +
               ((j % 2) * (num_slices[0] + 1)) + ((j / 2) * (num_slices[0] + 1) * (num_slices[1] + 1));
      }
    }

    /**
     * \brief Vertex at cell adjacencies
     *
     * \param[in] idx  Cell index
     * \param[in] vert Vertex index (0, 1, 2, 3, 4, 5, 6, 7)
     *
     * \returns The index of the j-th vertex of the i-th cell
     */
    std::size_t cell(std::size_t idx, std::size_t vert) const
    {
      const std::size_t z_layer = idx / ((_xs.size() - 1) * (_ys.size() - 1));
      const std::size_t y_layer = idx % ((_xs.size() - 1) * (_ys.size() - 1)) / (_xs.size() - 1);
      const std::size_t x_layer = idx % ((_xs.size() - 1) * (_ys.size() - 1)) % (_xs.size() - 1);

      const std::size_t base = (z_layer * (_xs.size() * _ys.size())) + (y_layer * _xs.size()) + x_layer;
      const std::size_t z_offset = _xs.size() * _ys.size();

      switch(vert)
      {
      case 0:
        return base;
      case 1:
        return base + 1;
      case 2:
        return base + _xs.size();
      case 3:
        return base + _xs.size() + 1;
      case 4:
        return base + z_offset;
      case 5:
        return base + 1 + z_offset;
      case 6:
        return base + _xs.size() + z_offset;
      case 7:
        return base + _xs.size() + 1 + z_offset;
      default:
        std::abort();
      }
    }

    /**
     * \brief Refinement steps for additional adaptive refinement
     *
     * \param[in] idx Vertex index
     *
     * \returns The recommended refinement level for the given vertex
     */
    std::uint64_t& subdivision_level(std::size_t idx)
    {
      return _subdivision_levels[idx];
    }

    /**
     * \brief Refinement steps for additional adaptive refinement
     *
     * \param[in] idx Vertex index
     *
     * \returns The recommended refinement level for the given vertex
     */
    const std::uint64_t& subdivision_level(std::size_t idx) const
    {
      return _subdivision_levels[idx];
    }

    /**
     * \brief Write mesh in the FEAT3 XML mesh file format
     *
     * \param[in] stream Outputstream to write to
     *
     * Writes this mesh in the FEAT3 XML mesh file format to the given output stream.
     * No refinements are applied to the mesh prior to writing.
     */
    void write_feat_xml(std::ostream& stream) const
    {
      stream << "<FeatMeshFile version=\"1\" mesh=\"conformal:hypercube:3:3\">\n";
      stream << "<Mesh type=\"conformal:hypercube:3:3\" size=\""
      << num_vertices() << " "
      << num_edges() << " "
      << num_faces() << " "
      << num_cells() << "\">\n";
      stream << "<Vertices>\n";
      for(std::size_t i(0); i < num_vertices(); i++)
      {
        const Point point = vertex(i);
        stream << point.x << " " << point.y << " " << point.z << "\n";
      }
      stream << "</Vertices>\n";

      stream << "<Topology dim=\"1\">\n";
      for(std::size_t i(0); i < num_edges(); i++)
      {
        for(std::size_t j(0); j < 2; j++)
        {
          stream << edge(i, j);
          stream << (j == 1 ? "\n" : " ");
        }
      }
      stream << "</Topology>\n";

      stream << "<Topology dim=\"2\">\n";
      for(std::size_t i(0); i < num_faces(); i++)
      {
        for(std::size_t j(0); j < 4; j++)
        {
          stream << face(i, j);
          stream << (j == 3 ? "\n" : " ");
        }
      }
      stream << "</Topology>\n";

      stream << "<Topology dim=\"3\">\n";
      for(std::size_t i(0); i < num_cells(); i++)
      {
        for(std::size_t j(0); j < 8; j++)
        {
          stream << cell(i, j);
          stream << (j == 7 ? "\n" : " ");
        }
      }
      stream << "</Topology>\n";

      stream << "</Mesh>\n";
      stream << "</FeatMeshFile>\n";
    }
  };

  /// Warning about intersecting triangles
  struct SelfIntersectionWarning
  {
    /// Name of this warning
    static constexpr std::string_view name = "self-intersection";

    /// First involved triangle
    std::uint32_t tri_a;

    /// Second involved triangle
    std::uint32_t tri_b;

    /// Constructor
    SelfIntersectionWarning(std::uint32_t a, std::uint32_t b) : tri_a(a), tri_b(b)
    {
    }
  };

  /// Warning about triangle with zero area
  struct DegenerateTriangleWarning
  {
    /// Name of this warning
    static constexpr std::string_view name = "degenerate-triangle";

    /// Index of degenerate triangle
    std::uint32_t idx;

    /// Constructor
    explicit DegenerateTriangleWarning(std::uint32_t i) : idx(i)
    {
    }
  };

  /// Warning about triangle with large aspect ratio
  struct AnisotropicTriangleWarning
  {
    /// Name of this warning
    static constexpr std::string_view name = "anisotropic-triangle";

    /// Index of anisotropic triangle
    std::uint32_t idx;

    /// Constructor
    explicit AnisotropicTriangleWarning(std::uint32_t i) : idx(i)
    {
    }
  };

  /// Collection of mesh warnings
  struct MeshWarnings
  {
    /// Collection of self-intersection warnings
    std::vector<SelfIntersectionWarning> self_intersections;

    /// Collection of degenerate-triangle warnings
    std::vector<DegenerateTriangleWarning> degenerate_triangles;

    /// Collection of anisotropic-triangle warnings
    std::vector<AnisotropicTriangleWarning> anisotropic_triangles;
  };

  /**
   * \brief Tagged union for error-handling. Contains either a success value of type T or an error value of type E
   */
  template<typename T, typename E>
  class Result
  {
    /**
     * \brief Type-tag for constructing a success-result.
     */
    struct TagOK
    {
    };

    /**
     * \brief Type-tag for constructing an error-result.
     */
    struct TagErr
    {
    };

    /// Result content
    std::variant<T, E> content;


    /// Ccopy constructor for success-results
    Result(TagOK /*tag*/, const T& value) : content(std::in_place_index<0>, value)
    {
    }

    /// Move constructor for success-results
    Result(TagOK /*tag*/, T&& value) : content(std::in_place_index<0>, std::move(value))
    {
    }

    /// Copy constructor for error-results
    Result(TagErr /*tag*/, const E& value) : content(std::in_place_index<1>, value)
    {
    }

    /// Move constructor for error-results
    Result(TagErr /*tag*/, E&& value) :
      content(std::in_place_index<1>, std::move(value))
    {
    }

  public:

    // Rule of five
    Result() = delete;

    /// Deleted copy constructor
    Result(const Result&) = delete;

    /// Deleted copy assignment
    Result& operator=(const Result&) = delete;

    /// Move constructor
    Result(Result&&) = default;

    /// Move-assignment
    Result& operator=(Result&&) = default;

    /// Destructor
    ~Result() = default;

    // Static factory methods

    /**
     * \brief Create a success-result from the given value
     *
     * \param[in] value Value of success
     */
    static Result ok(const T& value)
    {
      return Result(TagOK{}, value);
    }

    /**
     * \brief Create a success-result from the given value
     *
     * \param[in] value Value of success
     */
    static Result ok(T&& value)
    {
      return Result(TagOK{}, std::move(value));
    }

    /**
     * \brief Create a error-result from the given value
     *
     * \param[in] value Value of error
     */
    static Result err(const E& value)
    {
      return Result(TagErr{}, value);
    }

    /**
     * \brief Create a error-result from the given value
     *
     * \param[in] value Value of error
     */
    static Result err(E&& value)
    {
      return Result(TagErr{}, std::move(value));
    }

    /// Returns true if this is a success-result
    bool is_ok() const
    {
      return content.index() == 0;
    }

    /// Returns true if this is an error-result
    bool is_err() const
    {
      return content.index() == 1;
    }

    /**
     * \brief Success value accessor
     *
     * \returns A copy of the stored value. Errors if this is not a success-result.
     */
    T ok_value() const
    {
      return std::get<0>(content);
    }

    /**
     * \brief Error value accessor
     *
     * \returns A copy of the stored value. Errors if this is not an error-result.
     */
    E err_value() const
    {
      return std::get<1>(content);
    }

    /**
     * \brief Success reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not a success-result.
     */
    T& ok_ref()
    {
      return std::get<0>(content);
    }

    /**
     * \brief Error reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not an error-result.
     */
    E& err_ref()
    {
      return std::get<1>(content);
    }

    /**
     * \brief Success reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not a success-result.
     */
    const T& ok_ref() const
    {
      return std::get<0>(content);
    }

    /**
     * \brief Error reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not an error-result.
     */
    const E& err_ref() const
    {
      return std::get<1>(content);
    }

    /**
     * \brief Success move accessor
     *
     * \returns A rvalue-reference to the stored value. Errors if this is not a success-result.
     */
    T&& take_ok() &&
    {
      return std::get<0>(std::move(content));
    };

    /**
     * \brief Error move accessor
     *
     * \returns A rvalue-reference to the stored value. Errors if this is not an error-result.
     */
    E&& take_err() &&
    {
      return std::get<1>(std::move(content));
    };
  };

  /**
   * \brief Tagged union for error-handling. Contains either a success value of type T or an error value of type E
   *
   * This is a specialization for results containing either nothing or an error.
   */
  template<typename E>
  class Result<void, E>
  {
    /// Optional error value
    std::optional<E> error;

  public:
    // Rule of five
    Result() : error(std::nullopt)
    {
    }

    /// Constructor
    explicit Result(const E& value) : error(value)
    {
    }

    /// Constructor
    explicit Result(E&& value) : error(std::move(value))
    {
    }

    /// Copy constructor
    Result(const Result&) = delete;

    /// Copy-assignment
    Result& operator=(const Result&) = delete;

    /// Move-constructor
    Result(Result&&) = default;

    /// Move-assignment
    Result& operator=(Result&&) = default;

    /// Destructor
    ~Result() = default;

    // Static factory methods

    /**
     * \brief Create a success-result
     */
    static Result ok()
    {
      return Result();
    }

    /**
     * \brief Create an error-result from the given value
     */
    static Result err(const E& value)
    {
      return Result(value);
    }

    /**
     * \brief Create an error-result by moving the given value
     */
    static Result err(E&& value)
    {
      return Result(std::move(value));
    }

    /// Returns true if this is a success-result
    bool is_ok() const
    {
      return !error.has_value();
    }

    /// Returns true if this is an error-result
    bool is_err() const
    {
      return error.has_value();
    }

    /**
     * \brief Error value accessor
     *
     * \returns A copy of the stored value. Errors if this is not an error-result.
     */
    E err_value() const
    {
      return error.value();
    }

    /**
     * \brief Error reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not an error-result.
     */
    E& err_ref()
    {
      return error.value();
    }

    /**
     * \brief Error reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not an error-result.
     */
    const E& err_ref() const
    {
      return error.value();
    }

    /**
     * \brief Error move accessor
     *
     * \returns A rvalue-reference to the stored value. Errors if this is not an error-result.
     */
    E&& take_err() &&
    {
      return std::move(error).value();
    };
  };
} // namespace MeshHexer
