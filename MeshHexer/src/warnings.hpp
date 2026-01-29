#pragma once

#include <cgal_types.hpp>
#include <meshhexer/types.hpp>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

namespace MeshHexer
{
  template<typename OutputIterator>
  void self_intersection_warnings(const Mesh& mesh, OutputIterator output)
  {
    std::vector<std::pair<FaceIndex, FaceIndex>> self_intersections;
    CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(self_intersections));

    for(auto& pair : self_intersections)
    {
      *output = SelfIntersectionWarning(pair.first, pair.second);
    }
  }

  template<typename OutputIterator>
  void degenerate_triangles_warnings(const Mesh& mesh, OutputIterator output)
  {
    std::array<Point3D, 3> points{};
    for(FaceIndex f : mesh.faces())
    {
      int idx = 0;
      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(f)))
      {
        points.at(idx++) = mesh.point(v);
      }

      if(idx == 3 && CGAL::collinear(points[0], points[1], points[2]))
      {
        *output = DegenerateTriangleWarning(f);
      }
    }
  }

  template<typename OutputIterator>
  void anisotropic_triangle_warnings(const Mesh& mesh, OutputIterator output)
  {
    for(FaceIndex f : mesh.faces())
    {
      if(CGAL::Polygon_mesh_processing::face_aspect_ratio(f, mesh) > 15.0)
      {
        *output = AnisotropicTriangleWarning(f);
      }
    }
  }

  inline void create_warnings(const Mesh& mesh, MeshWarnings& warnings)
  {
    self_intersection_warnings(mesh, std::back_inserter(warnings.self_intersections));
    degenerate_triangles_warnings(mesh, std::back_inserter(warnings.degenerate_triangles));
    anisotropic_triangle_warnings(mesh, std::back_inserter(warnings.anisotropic_triangles));
  }
} // namespace MeshHexer
