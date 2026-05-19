#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/number_utils.h>

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <limits>
#include <cmath>

namespace PMP = CGAL::Polygon_mesh_processing;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using SurfaceMesh = CGAL::Surface_mesh<Kernel::Point_3>;
using SphereTraits = CGAL::Min_sphere_of_spheres_d_traits_3<Kernel, double, CGAL::Tag_true>;
using MinSphere = CGAL::Min_sphere_of_spheres_d<SphereTraits>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh>;
using AABBTraits = CGAL::AABB_traits<Kernel, Primitive>;
using AABBTree = CGAL::AABB_tree<AABBTraits>;
using SideTester = CGAL::Side_of_triangle_mesh<SurfaceMesh, Kernel>;

struct MeshHandle {
    SurfaceMesh mesh;
    mutable std::unique_ptr<AABBTree> tree;
    mutable std::unique_ptr<SideTester> side_tester;
};

namespace {

bool face_is_degenerate(const SurfaceMesh& mesh, SurfaceMesh::Face_index face) {
    Kernel::Point_3 points[3];
    std::size_t count = 0;
    SurfaceMesh::Vertex_index vertices[3];

    for (const auto vertex : CGAL::vertices_around_face(mesh.halfedge(face), mesh)) {
        if (count >= 3) {
            return true;
        }
        vertices[count] = vertex;
        points[count] = mesh.point(vertex);
        ++count;
    }

    if (count != 3) {
        return true;
    }

    if (vertices[0] == vertices[1] || vertices[0] == vertices[2] || vertices[1] == vertices[2]) {
        return true;
    }

    for (const auto& point : points) {
        if (!std::isfinite(point.x()) || !std::isfinite(point.y()) || !std::isfinite(point.z())) {
            return true;
        }
    }

    return CGAL::collinear(points[0], points[1], points[2]);
}

bool validate_mesh_for_queries(const SurfaceMesh& mesh, std::string& reason) {
    if (CGAL::is_empty(mesh) || mesh.number_of_vertices() == 0 || mesh.number_of_faces() == 0) {
        reason = "mesh is empty";
        return false;
    }

    if (!CGAL::is_triangle_mesh(mesh)) {
        reason = "mesh is not purely triangular after triangulation";
        return false;
    }

    if (!CGAL::is_closed(mesh)) {
        reason = "mesh is not closed; inside/outside tests require a closed surface";
        return false;
    }

    std::size_t degenerate_face_count = 0;
    for (const auto face : mesh.faces()) {
        if (face_is_degenerate(mesh, face)) {
            ++degenerate_face_count;
        }
    }
    if (degenerate_face_count > 0) {
        reason = "mesh contains " + std::to_string(degenerate_face_count) + " degenerate triangle(s)";
        return false;
    }

    std::vector<std::pair<SurfaceMesh::Face_index, SurfaceMesh::Face_index>> self_intersections;
    PMP::self_intersections(mesh, std::back_inserter(self_intersections));
    if (!self_intersections.empty()) {
        reason = "mesh contains " + std::to_string(self_intersections.size()) +
                 " self-intersection pair(s); overlapping component interfaces are not supported";
        return false;
    }

    return true;
}

AABBTree* ensure_tree(const MeshHandle* handle) {
    if (handle == nullptr) {
        return nullptr;
    }
    if (!handle->tree) {
        auto face_range = faces(handle->mesh);
        handle->tree = std::make_unique<AABBTree>(face_range.first, face_range.second, handle->mesh);
        handle->tree->accelerate_distance_queries();
    }
    return handle->tree.get();
}

SideTester* ensure_side_tester(const MeshHandle* handle) {
    if (handle == nullptr) {
        return nullptr;
    }
    if (!handle->side_tester) {
        handle->side_tester = std::make_unique<SideTester>(handle->mesh);
    }
    return handle->side_tester.get();
}

}  // namespace

extern "C" {

MeshHandle* cgal_load_off(const char* path) {
    if (path == nullptr) {
        return nullptr;
    }

    auto handle = std::make_unique<MeshHandle>();
    std::ifstream input(path);
    if (!input || !(input >> handle->mesh)) {
        return nullptr;
    }

    // Ensure purely triangular faces for predictable connectivity.
    PMP::triangulate_faces(handle->mesh);

    std::string validation_error;
    if (!validate_mesh_for_queries(handle->mesh, validation_error)) {
        std::cerr << "CGAL rejected OFF mesh '" << path << "': " << validation_error << "\n";
        return nullptr;
    }

    return handle.release();
}

void cgal_free_mesh(MeshHandle* handle) {
    delete handle;
}

std::size_t cgal_get_vertex_count(const MeshHandle* handle) {
    if (handle == nullptr) {
        return 0;
    }
    return handle->mesh.number_of_vertices();
}

std::size_t cgal_get_triangle_count(const MeshHandle* handle) {
    if (handle == nullptr) {
        return 0;
    }
    return handle->mesh.number_of_faces();
}

int cgal_copy_vertices(const MeshHandle* handle, double* out_buffer, std::size_t buffer_len) {
    if (handle == nullptr || out_buffer == nullptr) {
        return 1;
    }

    const std::size_t vertex_count = handle->mesh.number_of_vertices();
    if (buffer_len < vertex_count * 3) {
        return 2;
    }

    std::size_t offset = 0;
    for (const auto vertex : handle->mesh.vertices()) {
        const auto& point = handle->mesh.point(vertex);
        out_buffer[offset++] = point.x();
        out_buffer[offset++] = point.y();
        out_buffer[offset++] = point.z();
    }

    return 0;
}

int cgal_copy_triangles(const MeshHandle* handle, int32_t* out_buffer, std::size_t buffer_len) {
    if (handle == nullptr || out_buffer == nullptr) {
        return 1;
    }

    const std::size_t face_count = handle->mesh.number_of_faces();
    if (buffer_len < face_count * 3) {
        return 2;
    }

    std::size_t offset = 0;
    for (const auto face : handle->mesh.faces()) {
        std::size_t local_vertex = 0;
        for (const auto vertex :
             CGAL::vertices_around_face(handle->mesh.halfedge(face), handle->mesh)) {
            if (local_vertex >= 3) {
                return 3;  // Should not happen because of triangulation.
            }
            // Convert to 1-based indexing for easier Fortran consumption.
            out_buffer[offset++] = static_cast<int32_t>(vertex.idx()) + 1;
            ++local_vertex;
        }
        if (local_vertex != 3) {
            return 4;  // Invalid triangle.
        }
    }

    return 0;
}

int cgal_minimal_bounding_sphere(const double* points, std::size_t point_count, double* center_out,
                                 double* radius_out) {
    if (points == nullptr || center_out == nullptr || radius_out == nullptr) {
        return 1;
    }
    if (point_count == 0) {
        return 2;
    }

    std::vector<SphereTraits::Sphere> spheres;
    spheres.reserve(point_count);
    for (std::size_t i = 0; i < point_count; ++i) {
        const double* coord = points + (i * 3);
        spheres.emplace_back(SphereTraits::Point(coord[0], coord[1], coord[2]), 0.0);
    }

    MinSphere min_sphere(spheres.begin(), spheres.end());
    if (min_sphere.is_empty() || !min_sphere.is_valid()) {
        return 3;
    }

    auto center_it = min_sphere.center_cartesian_begin();
    for (int i = 0; i < 3; ++i) {
        center_out[i] = CGAL::to_double(center_it[i]);
    }

    const auto radius = min_sphere.radius();
    if (radius < 0) {
        return 4;
    }
    *radius_out = CGAL::to_double(radius);
    return 0;
}

int cgal_sphere_intersects_mesh(const MeshHandle* handle, const double* center, double radius,
                                int32_t* intersects_out) {
    if (handle == nullptr || center == nullptr || intersects_out == nullptr) {
        return 1;
    }
    if (radius < 0.0) {
        return 2;
    }

    AABBTree* tree = ensure_tree(handle);
    if (tree == nullptr) {
        return 3;
    }

    const Kernel::Point_3 center_point(center[0], center[1], center[2]);
    const double squared_distance = CGAL::to_double(tree->squared_distance(center_point));
    const double radius_squared = radius * radius;

    *intersects_out = (squared_distance <= radius_squared) ? 1 : 0;
    return 0;
}

int cgal_point_inside_mesh(const MeshHandle* handle, const double* point, int32_t* inside_out) {
    if (handle == nullptr || point == nullptr || inside_out == nullptr) {
        return 1;
    }

    SideTester* tester = ensure_side_tester(handle);
    if (tester == nullptr) {
        return 2;
    }

    const Kernel::Point_3 query(point[0], point[1], point[2]);
    const auto side = (*tester)(query);
    *inside_out = (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY) ? 1 : 0;
    return 0;
}

int cgal_closest_point_index_to_mesh(const MeshHandle* handle, const double* points, std::size_t point_count,
                                     int32_t* index_out) {
    if (handle == nullptr || points == nullptr || index_out == nullptr) {
        return 1;
    }
    if (point_count == 0) {
        return 2;
    }

    AABBTree* tree = ensure_tree(handle);
    if (tree == nullptr) {
        return 3;
    }

    double best_distance = std::numeric_limits<double>::infinity();
    std::size_t best_index = 0;
    for (std::size_t i = 0; i < point_count; ++i) {
        const double* coord = points + (i * 3);
        const Kernel::Point_3 query(coord[0], coord[1], coord[2]);
        const double sq_dist = CGAL::to_double(tree->squared_distance(query));
        if (sq_dist < best_distance) {
            best_distance = sq_dist;
            best_index = i;
        }
    }

    *index_out = static_cast<int32_t>(best_index + 1);  // return 1-based index
    return 0;
}

int cgal_segment_intersects_mesh(const MeshHandle* handle, const double* segment_points, int32_t* intersects_out) {
    if (handle == nullptr || segment_points == nullptr || intersects_out == nullptr) {
        return 1;
    }

    AABBTree* tree = ensure_tree(handle);
    if (tree == nullptr) {
        return 2;
    }

    const Kernel::Point_3 p0(segment_points[0], segment_points[1], segment_points[2]);
    const Kernel::Point_3 p1(segment_points[3], segment_points[4], segment_points[5]);
    const Kernel::Segment_3 segment(p0, p1);

    *intersects_out = tree->do_intersect(segment) ? 1 : 0;
    return 0;
}

int cgal_triangle_intersects_mesh(const MeshHandle* handle, const double* triangle_points, int32_t* intersects_out) {
    if (handle == nullptr || triangle_points == nullptr || intersects_out == nullptr) {
        return 1;
    }

    AABBTree* tree = ensure_tree(handle);
    if (tree == nullptr) {
        return 2;
    }

    const Kernel::Point_3 p0(triangle_points[0], triangle_points[1], triangle_points[2]);
    const Kernel::Point_3 p1(triangle_points[3], triangle_points[4], triangle_points[5]);
    const Kernel::Point_3 p2(triangle_points[6], triangle_points[7], triangle_points[8]);
    const Kernel::Triangle_3 triangle(p0, p1, p2);

    *intersects_out = tree->do_intersect(triangle) ? 1 : 0;
    return 0;
}

int cgal_build_aabb_tree(MeshHandle* handle) {
    if (handle == nullptr) {
        return 1;
    }
    AABBTree* tree = ensure_tree(handle);
    if (tree == nullptr) {
        return 2;
    }
    return 0;
}

int cgal_signed_distance_to_mesh(const MeshHandle* handle, const double* point, double* distance_out) {
    if (handle == nullptr || point == nullptr || distance_out == nullptr) {
        return 1;
    }

    AABBTree* tree = ensure_tree(handle);
    SideTester* tester = ensure_side_tester(handle);
    if (tree == nullptr || tester == nullptr) {
        return 2;
    }

    const Kernel::Point_3 query(point[0], point[1], point[2]);
    const double sq_dist = CGAL::to_double(tree->squared_distance(query));
    const double dist = std::sqrt(std::max(0.0, sq_dist));
    const auto side = (*tester)(query);
    const bool is_inside = (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY);
    *distance_out = is_inside ? dist : -dist;
    return 0;
}

}  // extern "C"
