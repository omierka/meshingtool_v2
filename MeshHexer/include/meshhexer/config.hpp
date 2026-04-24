#pragma once

#include <cstddef>

#include <meshhexer/types.hpp>

namespace MeshHexer
{
  struct MinGapConfig
  {
    double neighbor_diameter_ratio_flag = 2.0;
    double neighbor_normal_angle_limit_deg = 60.0;
    double small_angle_limit_deg = 4.0;
    double relative_min_search_radius = 0.005;
    double dihedral_angle_threshold_rad = 0.3;
    double face_aspect_ratio_limit = 5.0;
    double edge_ratio_limit = 5.0;
    double min_diameter_fraction = 1.0e-4;
    double gap_score_threshold = 0.95;
    double score_percentile_fallback = 0.10;
    double monitor_histogram_base = 3.0;
    double monitor_histogram_eps = 1.0e-12;
    double histogram_target_fraction = 0.001;
    double histogram_min_fraction_threshold = 0.001;
    double histogram_max_fraction_threshold = 0.80;
    std::size_t max_histogram_span_bins = 3;
    std::size_t min_histogram_span_bins = 1;
    double coarse_mesh_scaling = 1.3;
  };

  const MinGapConfig& min_gap_config();
  const RoundGeometryAnalysisConfig& round_geometry_analysis_config();
} // namespace MeshHexer
