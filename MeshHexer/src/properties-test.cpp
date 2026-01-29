#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <properties.hpp>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <vector>

using namespace MeshHexer;

TEST_CASE( "Maximal inscribed spheres", "[properties]") {
  Mesh mesh;

  std::vector<Point3D> vertices(8);
  vertices[0] = Point3D(-5, -100, -100);
  vertices[1] = Point3D( 5, -100, -100);
  vertices[2] = Point3D(-5,  100, -100);
  vertices[3] = Point3D( 5,  100, -100);

  vertices[4] = Point3D(-5, -100,  100);
  vertices[5] = Point3D( 5, -100,  100);
  vertices[6] = Point3D(-5,  100,  100);
  vertices[7] = Point3D( 5,  100,  100);

  std::vector<std::vector<std::size_t>> face_polygons(12);

  face_polygons[0] = {0, 1, 2};
  face_polygons[1] = {2, 1, 3};
  face_polygons[2] = {4, 5, 6};
  face_polygons[3] = {6, 5, 7};
  face_polygons[4] = {0, 1, 4};
  face_polygons[5] = {4, 1, 5};
  face_polygons[6] = {2, 3, 6};
  face_polygons[7] = {6, 3, 7};
  face_polygons[8] = {0, 2, 4};
  face_polygons[9] = {2, 4, 6};
  face_polygons[10] = {1, 3, 5};
  face_polygons[11] = {5, 3, 7};

  REQUIRE(CGAL::Polygon_mesh_processing::orient_polygon_soup(vertices, face_polygons));
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(vertices, face_polygons, mesh);

  AABBTree aabb_tree(mesh.faces_begin(), mesh.faces_end(), mesh);

  compute_vertex_normals(mesh);
  maximal_inscribed_spheres(mesh, aabb_tree);

  auto maybe_mis = mesh.property_map<FaceIndex, double>("f:MIS_diameter");
  REQUIRE(maybe_mis.has_value());
  auto mis = maybe_mis.value();

  // Large sides should have diameter of 10
  REQUIRE_THAT(mis[FaceIndex(8)], Catch::Matchers::WithinAbs(10.0, 0.1));
  REQUIRE_THAT(mis[FaceIndex(9)], Catch::Matchers::WithinAbs(10.0, 0.1));
  REQUIRE_THAT(mis[FaceIndex(10)], Catch::Matchers::WithinAbs(10.0, 0.1));
  REQUIRE_THAT(mis[FaceIndex(11)], Catch::Matchers::WithinAbs(10.0, 0.1));
}
