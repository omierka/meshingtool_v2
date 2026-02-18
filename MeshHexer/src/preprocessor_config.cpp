#include <meshhexer/config.hpp>

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <sstream>
#include <string>

namespace MeshHexer
{
  namespace
  {
    std::filesystem::path locate_config_file()
    {
      if(const char* env = std::getenv("PREPROCESSOR_CONFIG"); env != nullptr && *env != '\0')
      {
        return std::filesystem::path(env);
      }
      return std::filesystem::current_path() / "preprocessor.cfg";
    }

    std::string strip_inline_comment(std::string value)
    {
      const auto comment_pos = value.find_first_of("#;");
      if(comment_pos != std::string::npos)
      {
        value = value.substr(0, comment_pos);
      }
      const auto is_space = [](unsigned char ch) { return std::isspace(ch) != 0; };
      auto begin = std::find_if_not(value.begin(), value.end(), is_space);
      auto end = std::find_if_not(value.rbegin(), value.rend(), is_space).base();
      if(begin >= end)
      {
        return {};
      }
      return std::string(begin, end);
    }

    std::optional<std::string> get_raw(const boost::property_tree::ptree& tree, const std::string& key)
    {
      if(auto node = tree.get_optional<std::string>(key))
      {
        return strip_inline_comment(*node);
      }
      return std::nullopt;
    }

    template<typename T>
    bool parse_value(const std::string& raw, T& out)
    {
      std::istringstream ss(raw);
      ss >> out;
      return static_cast<bool>(ss);
    }

    template<typename T>
    void load_value(const boost::property_tree::ptree& tree, const std::string& key, T& target)
    {
      if(auto raw = get_raw(tree, key))
      {
        T parsed{};
        if(parse_value(*raw, parsed))
        {
          target = parsed;
        }
      }
    }

    MinGapConfig load_config()
    {
      MinGapConfig cfg;
      const auto path = locate_config_file();
      try
      {
        boost::property_tree::ptree tree;
        boost::property_tree::ini_parser::read_ini(path.string(), tree);

        load_value(tree, "2DMeshAnalysis/MingapValidity.neighbor_diameter_ratio_flag", cfg.neighbor_diameter_ratio_flag);
        load_value(tree, "2DMeshAnalysis/MingapValidity.neighbor_normal_angle_limit_deg", cfg.neighbor_normal_angle_limit_deg);
        load_value(tree, "2DMeshAnalysis/MingapValidity.small_angle_limit_deg", cfg.small_angle_limit_deg);

        load_value(tree, "2DMeshAnalysis/MingapScore.relative_min_search_radius", cfg.relative_min_search_radius);
        load_value(tree, "2DMeshAnalysis/MingapScore.dihedral_angle_threshold_rad", cfg.dihedral_angle_threshold_rad);
        load_value(tree, "2DMeshAnalysis/MingapScore.face_aspect_ratio_limit", cfg.face_aspect_ratio_limit);
        load_value(tree, "2DMeshAnalysis/MingapScore.edge_ratio_limit", cfg.edge_ratio_limit);
        load_value(tree, "2DMeshAnalysis/MingapScore.min_diameter_fraction", cfg.min_diameter_fraction);
        load_value(tree, "2DMeshAnalysis/MingapScore.gap_score_threshold", cfg.gap_score_threshold);
        load_value(tree, "2DMeshAnalysis/MingapScore.score_percentile_fallback", cfg.score_percentile_fallback);

        load_value(tree, "2DMeshAnalysis/MingapDesign.monitor_histogram_base", cfg.monitor_histogram_base);
        load_value(tree, "2DMeshAnalysis/MingapDesign.monitor_histogram_eps", cfg.monitor_histogram_eps);
        load_value(tree, "2DMeshAnalysis/MingapDesign.histogram_target_fraction", cfg.histogram_target_fraction);
        load_value(tree, "2DMeshAnalysis/MingapDesign.histogram_min_fraction_threshold", cfg.histogram_min_fraction_threshold);
        load_value(tree, "2DMeshAnalysis/MingapDesign.histogram_max_fraction_threshold", cfg.histogram_max_fraction_threshold);

        {
          std::size_t parsed = cfg.max_histogram_span_bins;
          if(auto raw = get_raw(tree, "2DMeshAnalysis/MingapDesign.max_histogram_span_bins"))
          {
            std::size_t temp;
            if(parse_value(*raw, temp))
            {
              parsed = temp;
            }
          }
          cfg.max_histogram_span_bins = parsed;
        }
        {
          std::size_t parsed = cfg.min_histogram_span_bins;
          if(auto raw = get_raw(tree, "2DMeshAnalysis/MingapDesign.min_histogram_span_bins"))
          {
            std::size_t temp;
            if(parse_value(*raw, temp))
            {
              parsed = temp;
            }
          }
          cfg.min_histogram_span_bins = parsed;
        }

        load_value(tree, "2DMeshAnalysis/MingapDesign.coarse_mesh_scaling", cfg.coarse_mesh_scaling);
      }
      catch(const std::exception&)
      {
        // keep defaults
      }
      return cfg;
    }
  } // namespace

  const MinGapConfig& min_gap_config()
  {
    static const MinGapConfig config = load_config();
    return config;
  }
} // namespace MeshHexer
