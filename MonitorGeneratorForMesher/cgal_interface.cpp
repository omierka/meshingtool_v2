#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/Bbox_3.h>
#include <vector>
#include <algorithm>
#include <memory>
#include <cstddef>

namespace {

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Triangle = Kernel::Triangle_3;
using Storage = std::vector<Triangle>;
using Iterator = Storage::iterator;
using Primitive = CGAL::AABB_triangle_primitive_3<Kernel, Iterator>;
using Traits = CGAL::AABB_traits_3<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

struct TriangleTreeHandle {
  Storage triangles;
  std::unique_ptr<Tree> tree;

  explicit TriangleTreeHandle(Storage&& tris)
      : triangles(std::move(tris)),
        tree(std::make_unique<Tree>(triangles.begin(), triangles.end())) {
    tree->accelerate_distance_queries();
  }
};

inline bool index_valid(int idx, int limit) {
  return idx >= 0 && idx < limit;
}

inline Point fetch_point(const double* coords, int vertex_index) {
  const double* base = coords + static_cast<std::size_t>(vertex_index) * 3u;
  return Point(base[0], base[1], base[2]);
}

}  // namespace

extern "C" {

TriangleTreeHandle* cgal_create_triangle_tree(const double* points,
                                              int n_vertices,
                                              const int* triangles,
                                              int n_triangles) {
  if (!points || !triangles || n_vertices <= 0 || n_triangles <= 0) {
    return nullptr;
  }

  Storage storage;
  storage.reserve(n_triangles);

  for (int i = 0; i < n_triangles; ++i) {
    const int* tri = triangles + 3 * static_cast<std::size_t>(i);
    const int id0 = tri[0] - 1;
    const int id1 = tri[1] - 1;
    const int id2 = tri[2] - 1;
    if (!index_valid(id0, n_vertices) ||
        !index_valid(id1, n_vertices) ||
        !index_valid(id2, n_vertices)) {
      continue;
    }
    storage.emplace_back(fetch_point(points, id0),
                         fetch_point(points, id1),
                         fetch_point(points, id2));
  }

  if (storage.empty()) {
    return nullptr;
  }

  try {
    return new TriangleTreeHandle(std::move(storage));
  } catch (...) {
    return nullptr;
  }
}

void cgal_free_triangle_tree(TriangleTreeHandle* handle) {
  delete handle;
}

int cgal_count_triangles_in_bbox(TriangleTreeHandle* handle,
                                 const double* bbox_min,
                                 const double* bbox_max) {
  if (!handle || !bbox_min || !bbox_max) {
    return 0;
  }
  const CGAL::Bbox_3 bbox(bbox_min[0], bbox_min[1], bbox_min[2],
                          bbox_max[0], bbox_max[1], bbox_max[2]);
  std::vector<Tree::Primitive_id> hits;
  handle->tree->all_intersected_primitives(bbox, std::back_inserter(hits));
  return static_cast<int>(hits.size());
}

int cgal_collect_triangles_in_bbox(TriangleTreeHandle* handle,
                                   const double* bbox_min,
                                   const double* bbox_max,
                                   int* out_indices,
                                   int max_indices) {
  if (!handle || !bbox_min || !bbox_max || !out_indices || max_indices <= 0) {
    return 0;
  }
  const CGAL::Bbox_3 bbox(bbox_min[0], bbox_min[1], bbox_min[2],
                          bbox_max[0], bbox_max[1], bbox_max[2]);
  std::vector<Tree::Primitive_id> hits;
  handle->tree->all_intersected_primitives(bbox, std::back_inserter(hits));
  const int available = static_cast<int>(hits.size());
  const int to_copy = std::min(available, max_indices);
  for (int i = 0; i < to_copy; ++i) {
    const Iterator it = hits[static_cast<std::size_t>(i)];
    const int idx = static_cast<int>(std::distance(handle->triangles.begin(), it));
    out_indices[i] = idx + 1;
  }
  return to_copy;
}

}  // extern "C"
