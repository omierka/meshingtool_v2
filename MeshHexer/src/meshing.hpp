#pragma once

#include <cgal_types.hpp>
#include <meshhexer/types.hpp>

namespace MeshHexer
{
  struct CuttingPlane
  {
    Plane3D plane;

    Point3D origin;
    Vector3D x_axis;
    Vector3D y_axis;

    Point2D project(const Point3D&) const;
  };

  class CrossSectionSampler
  {
  public:
    virtual ~CrossSectionSampler() = default;
    virtual Point3D origin() const = 0;
    virtual int num_planes() const = 0;
    virtual CuttingPlane get_plane(int) const = 0;
    virtual Point2D project(Point3D p) const = 0;
    virtual CuttingPlane get_plane_through_vertex(Point3D p) const = 0;
  };

  class RadialCrossSectionSampler : public CrossSectionSampler
  {
    int _num_planes;

    Point3D _origin;

    Vector3D _up;
    Vector3D _u;
    Vector3D _v;

  public:
    RadialCrossSectionSampler(int num_planes, const Point3D& o, const Vector3D& normal, const Vector3D& up) :
      _num_planes(num_planes),
      _origin(o),
      _up(up),
      _u(normal),
      _v(CGAL::cross_product(normal, up))
    {
      if(CGAL::scalar_product(normal, up) > 1e-8)
      {
        throw std::invalid_argument("RadialCrossSectionSampler: normal and up must be orthogonal!");
      }

      _up = _up / CGAL::approximate_sqrt(_up.squared_length());
      _u = _u / CGAL::approximate_sqrt(_u.squared_length());
      _v = _v / CGAL::approximate_sqrt(_v.squared_length());
    }

    Point3D origin() const override;
    int num_planes() const override;
    CuttingPlane get_plane(int idx) const override;

    Point2D project(Point3D p) const override;
    CuttingPlane get_plane_through_vertex(Point3D p) const override;
  };

  class LineCrossSectionSampler : public CrossSectionSampler
  {
    int _num_planes;

    Point3D _start;
    Point3D _end;

    Vector3D _normal;
    Vector3D _x_axis;
    Vector3D _y_axis;

  public:
    LineCrossSectionSampler(
      int num_planes,
      const Point3D& start,
      const Point3D& end,
      const Vector3D& normal,
      const Vector3D& up) :
      _num_planes(num_planes),
      _start(start),
      _end(end),
      _normal(normal),
      _x_axis(CGAL::cross_product(normal, up)),
      _y_axis(up)
    {
      if(CGAL::scalar_product(normal, up) > 1e-8)
      {
        throw std::invalid_argument("LineCrossSectionSampler: normal and up must be orthogonal!");
      }

      _normal = _normal / CGAL::approximate_sqrt(_normal.squared_length());
      _x_axis = _x_axis / CGAL::approximate_sqrt(_x_axis.squared_length());
      _y_axis = _y_axis / CGAL::approximate_sqrt(_y_axis.squared_length());
    }

    Point3D origin() const override;
    int num_planes() const override;
    CuttingPlane get_plane(int idx) const override;

    Point2D project(Point3D p) const override;
    CuttingPlane get_plane_through_vertex(Point3D p) const override;
  };

  std::pair<Polygon2D, std::vector<std::size_t>> simplify_by_normal(
    const Polygon2D& polygon,
    const std::function<bool(const std::vector<Vector2D>&, const Vector2D&)>& continue_pred);

  std::vector<PolygonWithHoles2D> union_of_cross_sections(const Mesh& mesh, const CrossSectionSampler& sampler);

  VolumeMesh fbm_mesh(Mesh& mesh, const FBMMeshSettings& settings);
} // namespace MeshHexer
