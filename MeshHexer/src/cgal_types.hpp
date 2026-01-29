#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangle_3.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>

// #include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>

namespace MeshHexer
{
  // Exact types
  using ExactKernel = CGAL::Exact_predicates_exact_constructions_kernel;

  // Inexact types

  using Kernel = CGAL::Simple_cartesian<double>;

  // using PrincipalVertexCurvature = CGAL::Polygon_mesh_processing::Principal_curvatures_and_directions<Kernel>;

  using Real = Kernel::RT;

  using Point3D = Kernel::Point_3;
  using Vector3D = Kernel::Vector_3;
  using Plane3D = Kernel::Plane_3;
  using Sphere3D = Kernel::Sphere_3;
  using Triangle3D = Kernel::Triangle_3;
  using Ray3D = Kernel::Ray_3;
  using Segment3D = Kernel::Segment_3;
  using Polyline3D = std::vector<Point3D>;
  using Polylines3D = std::list<Polyline3D>;

  using Point2D = Kernel::Point_2;
  using Vector2D = Kernel::Vector_2;
  using Polyline2D = std::vector<Point2D>;
  using Polylines2D = std::list<Polyline2D>;
  using Polygon2D = CGAL::Polygon_2<Kernel>;
  using PolygonWithHoles2D = CGAL::Polygon_with_holes_2<Kernel>;

  // Mesh types
  using Mesh = CGAL::Surface_mesh<Point3D>;
  using VertexIndex = Mesh::Vertex_index;
  using EdgeIndex = Mesh::Edge_index;
  using FaceIndex = Mesh::Face_index;
  using HalfedgeIndex = Mesh::Halfedge_index;

  // AABBTree types
  using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
  using AABBTraits = CGAL::AABB_traits_3<Kernel, Primitive>;
  using AABBTree = CGAL::AABB_tree<AABBTraits>;

  // Intersection types
  using RayIntersection = AABBTree::Intersection_and_primitive_id<Ray3D>::Type;
} // namespace MeshHexer
