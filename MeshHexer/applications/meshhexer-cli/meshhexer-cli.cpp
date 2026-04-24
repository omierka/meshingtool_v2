#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <unistd.h>

#include <cgal_types.hpp>
#include <meshhexer/meshhexer.hpp>
#include <meshhexer/types.hpp>
#include <meshhexer/config.hpp>
#include <properties.hpp>
#include <warnings.hpp>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace MeshHexerCLI::Markdown
{
  static std::string h1(const std::string& heading)
  {
    return heading + "\n" + std::string(heading.size(), '=');
  }

  static std::string h2(const std::string& heading)
  {
    return heading + "\n" + std::string(heading.size(), '-');
  }

  static std::string li(const std::string& content)
  {
    return "- " + content;
  }
} // namespace MeshHexerCLI::Markdown

namespace MeshHexerCLI
{
  namespace
  {
    constexpr double SMALL_FACE_RELATIVE_AREA = 1e-20;

    struct ParsedRoundSetup
    {
      std::string geometry_type;
      bool has_geometry_type = false;
      std::string geometry_start_raw;
      std::string geometry_length_raw;
      bool has_geometry_start = false;
      bool has_geometry_length = false;
      std::string preprocessing_geometry_start_raw;
      std::string preprocessing_geometry_length_raw;
      bool has_preprocessing_geometry_start = false;
      bool has_preprocessing_geometry_length = false;
      std::string preprocessing_hexmesher;
      bool has_preprocessing_hexmesher = false;
      std::vector<MeshHexer::RoundGeometryInflow> inflows;
    };

    struct SetupRewriteResult
    {
      std::vector<std::string> lines;
      std::size_t shifted_centers = 0;
      std::size_t shifted_midpoints = 0;
      bool preprocessing_section_updated = false;
    };

    std::string classification_to_string(MeshHexer::RoundGeometryClassification classification);

    std::string trim_copy(const std::string& value)
    {
      const auto begin = value.find_first_not_of(" \t\r\n");
      if(begin == std::string::npos)
      {
        return "";
      }
      const auto end = value.find_last_not_of(" \t\r\n");
      return value.substr(begin, end - begin + 1);
    }

    std::string lowercase_copy(std::string value)
    {
      std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      return value;
    }

    std::string uppercase_copy(std::string value)
    {
      std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
      return value;
    }

    std::optional<MeshHexer::Point> parse_point3(const std::string& raw)
    {
      std::stringstream stream(raw);
      std::string token;
      std::vector<double> values;
      while(std::getline(stream, token, ','))
      {
        token = trim_copy(token);
        if(token.empty())
        {
          return std::nullopt;
        }
        try
        {
          values.push_back(std::stod(token));
        }
        catch(const std::exception&)
        {
          return std::nullopt;
        }
      }
      if(values.size() != 3)
      {
        return std::nullopt;
      }
      return MeshHexer::Point{values[0], values[1], values[2]};
    }

    std::optional<std::size_t> parse_inflow_section_index(const std::string& section_name)
    {
      const std::string prefix = "e3dprocessparameters/inflow_";
      const std::string section = lowercase_copy(section_name);
      if(section.rfind(prefix, 0) != 0)
      {
        return std::nullopt;
      }
      const std::string suffix = section.substr(prefix.size());
      if(suffix.empty())
      {
        return std::nullopt;
      }
      try
      {
        const std::size_t value = static_cast<std::size_t>(std::stoul(suffix));
        if(value == 0)
        {
          return std::nullopt;
        }
        return value - 1;
      }
      catch(const std::exception&)
      {
        return std::nullopt;
      }
    }

    MeshHexer::Result<ParsedRoundSetup, std::string> parse_round_setup_file(const std::filesystem::path& setup_path)
    {
      using ResultType = MeshHexer::Result<ParsedRoundSetup, std::string>;

      std::ifstream input(setup_path);
      if(!input)
      {
        return ResultType::err("Failed to open setup file " + setup_path.string());
      }

      ParsedRoundSetup result;
      std::string current_section;
      std::string line;
      while(std::getline(input, line))
      {
        line = trim_copy(line);
        if(line.empty() || line[0] == '#' || line[0] == ';' || line[0] == '!')
        {
          continue;
        }
        if(line.front() == '[' && line.back() == ']')
        {
          current_section = trim_copy(line.substr(1, line.size() - 2));
          continue;
        }

        const auto eq_pos = line.find('=');
        if(eq_pos == std::string::npos || eq_pos == 0 || eq_pos + 1 >= line.size())
        {
          continue;
        }

        const std::string key = lowercase_copy(trim_copy(line.substr(0, eq_pos)));
        const std::string value = trim_copy(line.substr(eq_pos + 1));

        if(lowercase_copy(current_section) == "e3dgeometrydata/machine" && key == "geometrytype")
        {
          result.geometry_type = uppercase_copy(value);
          result.has_geometry_type = !result.geometry_type.empty();
          continue;
        }
        if(lowercase_copy(current_section) == "e3dgeometrydata/machine" && key == "geometrystart")
        {
          result.geometry_start_raw = value;
          result.has_geometry_start = !result.geometry_start_raw.empty();
          continue;
        }
        if(lowercase_copy(current_section) == "e3dgeometrydata/machine" && key == "geometrylength")
        {
          result.geometry_length_raw = value;
          result.has_geometry_length = !result.geometry_length_raw.empty();
          continue;
        }
        if(lowercase_copy(current_section) == "e3dgeometrydata/preprocessing" && key == "geometrystart")
        {
          result.preprocessing_geometry_start_raw = value;
          result.has_preprocessing_geometry_start = !result.preprocessing_geometry_start_raw.empty();
          continue;
        }
        if(lowercase_copy(current_section) == "e3dgeometrydata/preprocessing" && key == "geometrylength")
        {
          result.preprocessing_geometry_length_raw = value;
          result.has_preprocessing_geometry_length = !result.preprocessing_geometry_length_raw.empty();
          continue;
        }
        if(lowercase_copy(current_section) == "e3dgeometrydata/preprocessing" && key == "hexmesher")
        {
          result.preprocessing_hexmesher = value;
          result.has_preprocessing_hexmesher = !result.preprocessing_hexmesher.empty();
          continue;
        }

        const auto maybe_index = parse_inflow_section_index(current_section);
        if(!maybe_index.has_value())
        {
          continue;
        }

        const std::size_t index = maybe_index.value();
        if(result.inflows.size() <= index)
        {
          result.inflows.resize(index + 1);
        }

        if(key == "center")
        {
          const auto point = parse_point3(value);
          if(point.has_value())
          {
            result.inflows[index].center = point.value();
            result.inflows[index].has_center = true;
          }
        }
        else if(key == "normal")
        {
          const auto point = parse_point3(value);
          if(point.has_value())
          {
            result.inflows[index].normal = point.value();
            result.inflows[index].has_normal = true;
          }
        }
      }

      result.inflows.erase(
        std::remove_if(
          result.inflows.begin(),
          result.inflows.end(),
          [](const MeshHexer::RoundGeometryInflow& inflow) { return !(inflow.has_center && inflow.has_normal); }),
        result.inflows.end());

      if(result.inflows.empty() && result.geometry_type == "ROUND")
      {
        return ResultType::err("No valid inflows with center and normal were found in " + setup_path.string());
      }

      return ResultType::ok(std::move(result));
    }

    MeshHexer::Result<std::filesystem::path, std::string> resolve_report_mesh_path(
      const std::filesystem::path& input_path)
    {
      using ResultType = MeshHexer::Result<std::filesystem::path, std::string>;

      std::error_code ec;
      const std::filesystem::file_status status = std::filesystem::status(input_path, ec);
      if(ec || !std::filesystem::exists(status))
      {
        return ResultType::err("Input path does not exist: " + input_path.string());
      }

      if(std::filesystem::is_directory(status))
      {
        const std::filesystem::path surface_path = input_path / "surface.off";
        const std::filesystem::file_status surface_status = std::filesystem::status(surface_path, ec);
        if(ec || !std::filesystem::exists(surface_status))
        {
          return ResultType::err("Expected surface.off in case directory " + input_path.string());
        }
        return ResultType::ok(surface_path);
      }

      return ResultType::err(
        "report expects a case directory containing surface.off and setup.e3d: " + input_path.string());
    }

    MeshHexer::Result<std::filesystem::path, std::string> resolve_case_directory(
      const std::filesystem::path& input_path)
    {
      using ResultType = MeshHexer::Result<std::filesystem::path, std::string>;

      std::error_code ec;
      const std::filesystem::file_status status = std::filesystem::status(input_path, ec);
      if(ec || !std::filesystem::exists(status) || !std::filesystem::is_directory(status))
      {
        return ResultType::err("Expected a case directory: " + input_path.string());
      }

      const std::filesystem::path setup_path = input_path / "setup.e3d";
      const std::filesystem::file_status setup_status = std::filesystem::status(setup_path, ec);
      if(ec || !std::filesystem::exists(setup_status))
      {
        return ResultType::err("Expected setup.e3d in case directory " + input_path.string());
      }

      return ResultType::ok(input_path);
    }

    std::string format_point3_xy_rounded(const MeshHexer::Point& point, std::optional<int> xy_decimals = std::nullopt)
    {
      std::ostringstream out;
      if(xy_decimals.has_value())
      {
        out << std::fixed << std::setprecision(*xy_decimals) << point.x << "," << point.y;
      }
      else
      {
        out << std::setprecision(17) << point.x << "," << point.y;
      }
      out << "," << point.z;
      return out.str();
    }

    std::string format_scalar(double value, int decimals)
    {
      std::ostringstream out;
      out << std::fixed << std::setprecision(decimals) << value;
      return out.str();
    }

    std::vector<std::string> preprocessing_section_lines(
      const MeshHexer::RoundGeometryAnalysisResult& analysis,
      const MeshHexer::RoundGeometryAnalysisConfig& config)
    {
      const int decimals = std::max(0, config.preprocessing_export_precision);
      std::vector<std::string> lines;
      lines.push_back("BarrelDiameter = " + format_scalar(analysis.outer_diameter, decimals));
      if(analysis.classification == MeshHexer::RoundGeometryClassification::HollowCylinder)
      {
        lines.push_back("InnerDiameter = " + format_scalar(analysis.inner_diameter, decimals));
      }
      lines.push_back("BarrelLength = " + format_scalar(analysis.z_max_physical - analysis.z_min_physical, decimals));
      lines.push_back("AxialStartPosition = " + format_scalar(analysis.z_min_physical, decimals));
      if(analysis.classification == MeshHexer::RoundGeometryClassification::FullCylinder)
      {
        lines.push_back("");
        lines.push_back("HexMesher=" + classification_to_string(analysis.classification));
        lines.push_back("FullCylinderPeriodicity = " + std::to_string(config.preprocessing_fullcylinder_periodicity));
      }
      else
      {
        lines.push_back("");
        lines.push_back("HexMesher=" + classification_to_string(analysis.classification));
      }
      lines.push_back("sEl_Tangential = " + format_scalar(config.preprocessing_sel_tangential, decimals));
      lines.push_back("sEl_Radial = " + format_scalar(config.preprocessing_sel_radial, decimals));
      lines.push_back("sEl_Axial = " + format_scalar(config.preprocessing_sel_axial, decimals));
      return lines;
    }

    std::vector<std::string> box_preprocessing_section_lines(
      const ParsedRoundSetup& setup,
      const MeshHexer::RoundGeometryAnalysisConfig& config)
    {
      std::vector<std::string> lines;
      lines.push_back("geometryStart = " + setup.geometry_start_raw);
      lines.push_back("geometryLength = " + setup.geometry_length_raw);
      lines.push_back("");
      lines.push_back("HexMesher=Box");
      lines.push_back("sEl_x = " + format_scalar(config.preprocessing_sel_x, 1));
      lines.push_back("sEl_y = " + format_scalar(config.preprocessing_sel_y, 1));
      lines.push_back("sEl_z = " + format_scalar(config.preprocessing_sel_z, 1));
      return lines;
    }

    bool has_preconfigured_box_preprocessing(const ParsedRoundSetup& setup)
    {
      return setup.has_preprocessing_geometry_start &&
             setup.has_preprocessing_geometry_length &&
             setup.has_preprocessing_hexmesher &&
             lowercase_copy(trim_copy(setup.preprocessing_hexmesher)) == "box";
    }

    std::string split_inline_comment(std::string& value)
    {
      const std::size_t comment_pos = value.find_first_of("#;");
      if(comment_pos == std::string::npos)
      {
        return {};
      }
      std::string comment = value.substr(comment_pos);
      value = value.substr(0, comment_pos);
      return comment;
    }

    std::unordered_map<std::string, std::string> parse_section_key_values(
      const std::filesystem::path& setup_path,
      const std::string& target_section)
    {
      std::unordered_map<std::string, std::string> values;
      std::ifstream input(setup_path);
      if(!input)
      {
        return values;
      }

      std::string current_section;
      std::string line;
      while(std::getline(input, line))
      {
        const std::string trimmed_line = trim_copy(line);
        if(trimmed_line.empty() || trimmed_line[0] == '#' || trimmed_line[0] == ';' || trimmed_line[0] == '!')
        {
          continue;
        }
        if(trimmed_line.front() == '[' && trimmed_line.back() == ']')
        {
          current_section = lowercase_copy(trim_copy(trimmed_line.substr(1, trimmed_line.size() - 2)));
          continue;
        }
        if(current_section != lowercase_copy(target_section))
        {
          continue;
        }

        const std::size_t eq_pos = line.find('=');
        if(eq_pos == std::string::npos || eq_pos == 0 || eq_pos + 1 >= line.size())
        {
          continue;
        }

        std::string value = line.substr(eq_pos + 1);
        split_inline_comment(value);
        values[lowercase_copy(trim_copy(line.substr(0, eq_pos)))] = trim_copy(value);
      }

      return values;
    }

    bool round_preprocessing_section_matches(
      const std::filesystem::path& setup_path,
      const MeshHexer::RoundGeometryAnalysisResult& analysis,
      const MeshHexer::RoundGeometryAnalysisConfig& config)
    {
      static const std::vector<std::string> managed_keys = {
        "barreldiameter",
        "innerdiameter",
        "barrellength",
        "axialstartposition",
        "hexmesher",
        "fullcylinderperiodicity",
        "sel_tangential",
        "sel_radial",
        "sel_axial",
      };

      const std::unordered_map<std::string, std::string> existing =
        parse_section_key_values(setup_path, "E3DGeometryData/Preprocessing");

      std::unordered_map<std::string, std::string> expected;
      for(const std::string& line : preprocessing_section_lines(analysis, config))
      {
        const std::size_t eq_pos = line.find('=');
        if(eq_pos == std::string::npos || eq_pos == 0 || eq_pos + 1 >= line.size())
        {
          continue;
        }
        expected[lowercase_copy(trim_copy(line.substr(0, eq_pos)))] = trim_copy(line.substr(eq_pos + 1));
      }

      for(const std::string& key : managed_keys)
      {
        const auto expected_it = expected.find(key);
        const auto existing_it = existing.find(key);
        if(expected_it == expected.end())
        {
          if(existing_it != existing.end())
          {
            return false;
          }
          continue;
        }
        if(existing_it == existing.end() || existing_it->second != expected_it->second)
        {
          return false;
        }
      }

      return true;
    }

    MeshHexer::Result<SetupRewriteResult, std::string> rewrite_setup_for_round_geometry(
      const std::filesystem::path& setup_path,
      double shift_x,
      double shift_y,
      const MeshHexer::RoundGeometryAnalysisResult& analysis,
      const MeshHexer::RoundGeometryAnalysisConfig& config)
    {
      using ResultType = MeshHexer::Result<SetupRewriteResult, std::string>;

      std::ifstream input(setup_path);
      if(!input)
      {
        return ResultType::err("Failed to open setup file " + setup_path.string());
      }

      SetupRewriteResult result;
      std::string current_section;
      const std::string preprocessing_section = "e3dgeometrydata/preprocessing";
      const std::vector<std::string> preprocessing_lines = preprocessing_section_lines(analysis, config);
      bool preprocessing_section_found = false;
      auto append_preprocessing_lines = [&]() {
        result.lines.insert(result.lines.end(), preprocessing_lines.begin(), preprocessing_lines.end());
        result.preprocessing_section_updated = true;
      };
      std::string line;
      while(std::getline(input, line))
      {
        std::string updated_line = line;
        const std::string trimmed_line = trim_copy(line);

        if(!trimmed_line.empty() && trimmed_line.front() == '[' && trimmed_line.back() == ']')
        {
          if(lowercase_copy(current_section) == preprocessing_section)
          {
            append_preprocessing_lines();
          }
          current_section = trim_copy(trimmed_line.substr(1, trimmed_line.size() - 2));
          if(lowercase_copy(current_section) == preprocessing_section)
          {
            preprocessing_section_found = true;
          }
          result.lines.push_back(std::move(updated_line));
          continue;
        }

        const std::string current_section_lower = lowercase_copy(current_section);
        if(current_section_lower == preprocessing_section)
        {
          const std::size_t eq_pos = line.find('=');
          if(eq_pos != std::string::npos && eq_pos > 0 && eq_pos + 1 < line.size())
          {
            const std::string key = lowercase_copy(trim_copy(line.substr(0, eq_pos)));
            if(key == "barreldiameter" || key == "innerdiameter" || key == "barrellength" ||
               key == "axialstartposition" || key == "fullcylinderperiodicity" || key == "hexmesher" || key == "sel_tangential" ||
               key == "sel_radial" || key == "sel_axial")
            {
              result.preprocessing_section_updated = true;
              continue;
            }
          }
        }
        else if(parse_inflow_section_index(current_section).has_value())
        {
          const std::size_t eq_pos = line.find('=');
          if(eq_pos != std::string::npos && eq_pos > 0 && eq_pos + 1 < line.size())
          {
            const std::string key = lowercase_copy(trim_copy(line.substr(0, eq_pos)));
            if(key == "center" || key == "midpointa" || key == "midpointb")
            {
              std::string value = line.substr(eq_pos + 1);
              const std::string comment = split_inline_comment(value);
              const auto point = parse_point3(value);
              if(!point.has_value())
              {
                return ResultType::err(
                  "Failed to parse " + key + " in inflow section [" + current_section + "] of " + setup_path.string());
              }

              MeshHexer::Point shifted = point.value();
              shifted.x += shift_x;
              shifted.y += shift_y;

              updated_line = line.substr(0, eq_pos + 1) + " " + format_point3_xy_rounded(shifted, 2);
              if(!comment.empty())
              {
                updated_line += " " + comment;
              }

              if(key == "center")
              {
                ++result.shifted_centers;
              }
              else
              {
                ++result.shifted_midpoints;
              }
            }
          }
        }

        result.lines.push_back(std::move(updated_line));
      }

      if(lowercase_copy(current_section) == preprocessing_section)
      {
        append_preprocessing_lines();
      }
      else if(!preprocessing_section_found)
      {
        if(!result.lines.empty() && !result.lines.back().empty())
        {
          result.lines.push_back("");
        }
        result.lines.push_back("[E3DGeometryData/Preprocessing]");
        append_preprocessing_lines();
      }

      return ResultType::ok(std::move(result));
    }

    MeshHexer::Result<SetupRewriteResult, std::string> rewrite_setup_for_box_geometry(
      const std::filesystem::path& setup_path,
      const ParsedRoundSetup& setup,
      const MeshHexer::RoundGeometryAnalysisConfig& config)
    {
      using ResultType = MeshHexer::Result<SetupRewriteResult, std::string>;

      std::ifstream input(setup_path);
      if(!input)
      {
        return ResultType::err("Failed to open setup file " + setup_path.string());
      }
      if(!setup.has_geometry_start || !setup.has_geometry_length)
      {
        return ResultType::err(
          "BOX configuration requires geometryStart and geometryLength in [E3DGeometryData/Machine] of " +
          setup_path.string());
      }

      SetupRewriteResult result;
      std::string current_section;
      const std::string machine_section = "e3dgeometrydata/machine";
      const std::string preprocessing_section = "e3dgeometrydata/preprocessing";
      const std::vector<std::string> preprocessing_lines = box_preprocessing_section_lines(setup, config);
      bool preprocessing_section_found = false;
      auto append_preprocessing_lines = [&]() {
        result.lines.insert(result.lines.end(), preprocessing_lines.begin(), preprocessing_lines.end());
        result.preprocessing_section_updated = true;
      };

      std::string line;
      while(std::getline(input, line))
      {
        std::string updated_line = line;
        const std::string trimmed_line = trim_copy(line);

        if(!trimmed_line.empty() && trimmed_line.front() == '[' && trimmed_line.back() == ']')
        {
          if(lowercase_copy(current_section) == preprocessing_section)
          {
            append_preprocessing_lines();
          }
          current_section = trim_copy(trimmed_line.substr(1, trimmed_line.size() - 2));
          if(lowercase_copy(current_section) == preprocessing_section)
          {
            preprocessing_section_found = true;
          }
          result.lines.push_back(std::move(updated_line));
          continue;
        }

        const std::string current_section_lower = lowercase_copy(current_section);
        const std::size_t eq_pos = line.find('=');
        if(eq_pos != std::string::npos && eq_pos > 0 && eq_pos + 1 < line.size())
        {
          const std::string key = lowercase_copy(trim_copy(line.substr(0, eq_pos)));
          if(current_section_lower == machine_section &&
             (key == "geometrystart" || key == "geometrylength"))
          {
            continue;
          }
          if(current_section_lower == preprocessing_section &&
             (key == "geometrystart" || key == "geometrylength" || key == "hexmesher" ||
              key == "sel_x" || key == "sel_y" || key == "sel_z"))
          {
            result.preprocessing_section_updated = true;
            continue;
          }
        }

        result.lines.push_back(std::move(updated_line));
      }

      if(lowercase_copy(current_section) == preprocessing_section)
      {
        append_preprocessing_lines();
      }
      else if(!preprocessing_section_found)
      {
        if(!result.lines.empty() && !result.lines.back().empty())
        {
          result.lines.push_back("");
        }
        result.lines.push_back("[E3DGeometryData/Preprocessing]");
        append_preprocessing_lines();
      }

      return ResultType::ok(std::move(result));
    }

    bool write_text_lines(const std::filesystem::path& path, const std::vector<std::string>& lines)
    {
      std::ofstream output(path);
      if(!output)
      {
        return false;
      }

      for(const std::string& line : lines)
      {
        output << line << "\n";
      }
      return output.good();
    }

    std::string classification_to_string(MeshHexer::RoundGeometryClassification classification)
    {
      switch(classification)
      {
      case MeshHexer::RoundGeometryClassification::FullCylinder: return "FullCylinder";
      case MeshHexer::RoundGeometryClassification::HollowCylinder: return "HollowCylinder";
      case MeshHexer::RoundGeometryClassification::Unknown:
      default: return "Unknown";
      }
    }

    void print_min_gap(const MeshHexer::Gap& min_gap, bool verbose, double coarse_mesh_size)
    {
      if(verbose)
      {
        std::cout << "Min-gap of " << min_gap.diameter << " between faces " << min_gap.face << " and " << min_gap.opposite_face << "\n";
        std::cout << "Histogram span (bins): " << min_gap.bin_span << "\n";
        std::cout << "Suggested coarse mesh size: " << coarse_mesh_size << "\n";
        std::cout << "Use `SelectIDs(IDs=[0, " << min_gap.face << ", 0, " << min_gap.opposite_face
              << "], FieldType='CELL')` to select the chosen triangles in ParaView\n";
      }
      else
      {
        std::cout << min_gap.diameter << " " << min_gap.bin_span << " " << coarse_mesh_size << "\n";
      }
    }

    template<typename Iter>
    void write_range_as_mtx(std::ostream& stream, Iter begin, Iter end)
    {
      // NOTE(mmuegge): For compatability with FEAT3 we write everything as
      // real. We could inspect the type produced by the iterator and set
      // integer as type if appropriate, but then we would also need to update
      // the parsing logic in FEAT3.
      stream << "%%MatrixMarket matrix array real general\n";

      const auto size = std::distance(begin, end);
      stream << size << " 1\n";

      for(Iter it = begin; it != end; it++)
      {
        stream << *it << "\n";
      }
    }
  } // namespace


  //////////////////
  // Usage strings
  //////////////////

  const static char* const usage =
    "meshhexer-cli: Commandline tool for the MeshHexer library\n"
    "\n"
    "Usage:\n"
    "meshhexer-cli [<global args>] <command>\n"
    "\n"
    "Global options:\n"
    "\t-h, --help\n"
    "\t\tProduce this help text\n"
    "\t--checkpoint-path\n"
    "\t\tWrite checkpoint file with intermediary values. File extension must be .ply or .vtu\n"
    "\n"
    "Commands:\n"
    "\tfbm-mesh\n"
    "\t\tGenerate a non-fitting volume mesh from the surface mesh.\n"
    "\t\tThe mesh is output as fbm_mesh.xml in FEAT3's mesh file format.\n"
    "\t\tA separate fbm_mesh.mtx file is created with recommended adaptive\n"
    "\t\trefinement levels for all vertices.\n"
    "\t\tThe output mesh is constructed such that all cells are about the same\n"
    "\t\tsize as their local min-gaps.\n"
    "\n"
"\tmin-gap\n"
"\t\tCalculate smallest inside gap between opposite faces of the mesh\n"
"\n"
"\treport\n"
"\t\tPrint information about the mesh\n"
"\n"
"\twarnings\n"
"\t\tPrint warnings about the mesh. Warns about self-intersections,\n"
"\t\tdegenerate triangles, and anisotropic triangles.\n"
"\n"
"See meshhexer-cli <command> --help for more details on the commands.\n";

  const static char* const mingap_usage =
    "Usage: meshhexer-cli min-gap [<args>] <mesh>\n"
    "\n"
    "Compute the smallest inside gap between opposite faces of the mesh.\n"
    "\n"
    "Options:\n"
    "\t-h, --help\n"
    "\t\tProduce this help text\n"
    "\t--verbose\n"
    "\t\tPrint involved faces and score alongside the min-gap\n";

  const static char* const fbm_usage =
    "Usage: meshhexer-cli fbm-mesh [<args>] <mesh>\n"
    "\n"
    "Generate a non-fitting volume mesh from the surface mesh.\n"
    "The mesh is output in FEAT3's mesh file format."
    "A separate .mtx file is created with recommended adaptive\n"
    "refinement levels for all vertices.\n"
    "The output mesh is constructed such that all cells are about the same\n"
    "size as their local min-gaps.\n"
    "\n"
    "Options:\n"
    "\t-h, --help\n"
    "\t\tProduce this help text\n"
    "\t--levels\n"
    "\t\tSet size of multigrid-hierarchy. If passed, the mesh will be constructed\n"
    "\t\tsuch that the finest level of the hierarchy matches the min-gaps.\n"
    "\t--bounding-box\n"
    "\t\tSet a custom bounding box for the fbm mesh. If not set the bounding box of the surface mesh is used.\n"
    "\t--output\n"
    "\t\tSet filename for output files. Default is fbm_mesh.\n";

  const static char* const report_usage =
    "Usage: meshhexer-cli report [<args>] <case-dir>\n"
    "\n"
    "Print information about the case mesh\n"
    "\n"
    "Options:\n"
    "\t--configure-case-for-preprocessing\n"
    "\t\tRead GeometryType from setup.e3d and run the matching preprocessing configuration path.\n"
    "\t\tCurrently ROUND and BOX are implemented\n";

  const static char* const warnings_usage =
    "Usage: meshhexer-cli warnings [<args>] <mesh>\n"
    "\n"
    "Print warnings about the mesh. Warns about self-intersections,\n"
    "degenerate triangles, and anisotropic triangles.\n"
    "\n"
    "Options:\n"
    "\t--summarize\n"
    "\t\tSummarize warnings\n";

  /////////////////////
  // Parameter structs
  /////////////////////

  struct GlobalParameters
  {
    /// If true, show help text and end program
    bool show_help = false;

    /// Checkpoint file location. Checkpoint file is written if path is not empty.
    std::filesystem::path checkpoint_path;

    /// Command to run
    std::string command;
  };

  struct MinGapParameters
  {
    /// If true, show help text and end program
    bool show_help = false;

    /// If true, output involved faces, score, and paraview script along with min-gap
    bool verbose = false;

    /// Mesh file to calculate mingap of
    std::filesystem::path mesh_file;
  };

  struct FbmMeshParameters
  {
    /// If true, show help text and end program
    bool show_help = false;

    /// Number of levels of intended multigrid-hierarchy
    std::size_t levels = 0;

    /// Mesh file to create base mesh for
    std::filesystem::path mesh_file;

    /// Filename for output files
    std::string output = "fbm_mesh";

    /// Custom bounding box
    std::optional<MeshHexer::BoundingBox> bounding_box;
  };

  struct ReportParameters
  {
    /// If true, show help text and end program
    bool show_help = false;

    /// If true, configure the case for preprocessing based on GeometryType in setup.e3d
    bool configure_case_for_preprocessing = false;

    /// Mesh file path
    std::filesystem::path mesh_file;
  };

  struct WarningsParameters
  {
    /// If true, show help text and end program
    bool show_help = false;

    /// Summarize warnings
    bool summarize = false;

    /// Mesh file path
    std::filesystem::path mesh_file;
  };

  ////////////////
  // Arg parsing
  ////////////////

  static bool cmp_argument(const char* parameter, char* arg)
  {
    std::size_t n = std::min(std::strlen(parameter), std::strlen(arg));
    return std::strncmp(parameter, arg, n) == 0;
  }

  /// Decrement argc and move argv to next argument
  static void consume_arg(int& argc, char*** argv)
  {
    argc--;
    (*argv)++;
  }

  /// Extract argument. Consumes args as required.
  static char* parse_argument(const char* parameter, int& argc, char*** argv)
  {
    char* arg = (*argv)[0];

    while((*arg != 0) && *arg != '=')
    {
      arg++;
    }

    if(*arg == 0)
    {
      // --param arg case

      if(argc < 2)
      {
        std::cerr << "Error parsing parameter " << parameter << ". Expected additional argument\n";
        std::cerr << usage;
        std::exit(1);
      }

      consume_arg(argc, argv);
      char* result = (*argv)[0];
      consume_arg(argc, argv);

      return result;
    }
    else
    {
      // --param=arg case

      // We used up one argument
      consume_arg(argc, argv);

      return arg + 1;
    }
  }

  /// Extract arguments. Consumes args as required.
  static std::vector<char*> parse_arguments(const char* parameter, int& argc, char*** argv, int length)
  {
    char* arg = (*argv)[0];

    while((*arg != 0) && *arg != '=')
    {
      arg++;
    }

    if(*arg != 0)
    {
      std::cerr << "Error parsing parameter " << parameter << ". <param>=<args> syntax is not supported for parameters expecting multiple arguments\n";
      std::cerr << usage;
      std::exit(1);
    }
      // --param arg case

    if(argc < length + 1)
    {
      std::cerr << "Error parsing parameter " << parameter << ". Not enough arguments given\n";
      std::cerr << usage;
      std::exit(1);
    }

    consume_arg(argc, argv);

    std::vector<char*> result;

    for(int i(0); i < length; i++)
    {
      char* next = (*argv)[0];

      if(cmp_argument("--", next))
      {
        std::cerr << "Error parsing parameter " << parameter << ". Expected argument but found next parameter " << next << "\n";
        std::cerr << usage;
        std::exit(1);
      }

      result.push_back(next);
      consume_arg(argc, argv);
    }

    return result;
  }

  /**
   * \brief Parse global parameters
   */
  static MeshHexer::Result<GlobalParameters, std::string> parse_global_args(int& argc, char*** argv)
  {
    using Result = MeshHexer::Result<GlobalParameters, std::string>;

    GlobalParameters result;

    while(argc > 0)
    {
      char* cmd = (*argv)[0];

      if(cmp_argument("--help", cmd) || cmp_argument("-h", cmd))
      {
        consume_arg(argc, argv);
        result.show_help = true;
      }
      else if(cmp_argument("--checkpoint-path", cmd))
      {
        char* path = parse_argument("--checkpoint-path", argc, argv);
        result.checkpoint_path = path;
      }
      else
      {
        // Unknown argument
        break;
      }
    }

    if(argc > 0)
    {
      result.command = (*argv[0]);
      consume_arg(argc, argv);

      if(
        result.command != "fbm-mesh" && result.command != "min-gap" && result.command != "report" &&
        result.command != "warnings")
      {
        return Result::err("Invalid command " + result.command);
      }
    }
    else if(!result.show_help)
    {
      return Result::err("Expected command!");
    }

    return Result::ok(result);
  }

  /**
   * \brief Parse mingap parameters
   */
  static MeshHexer::Result<MinGapParameters, std::string> parse_mingap_args(int& argc, char*** argv)
  {
    using Result = MeshHexer::Result<MinGapParameters, std::string>;

    MinGapParameters result;

    while(argc > 0)
    {
      char* cmd = (*argv)[0];

      if(cmp_argument("--help", cmd) || cmp_argument("-h", cmd))
      {
        consume_arg(argc, argv);
        result.show_help = true;
      }
      else if(cmp_argument("--verbose", cmd))
      {
        consume_arg(argc, argv);
        result.verbose = true;
      }
      else
      {
        // Unknown argument
        break;
      }
    }

    if(argc > 0)
    {
      result.mesh_file = (*argv[0]);
      consume_arg(argc, argv);
    }
    else if(!result.show_help)
    {
      return Result::err("Expected mesh file!");
    }

    return Result::ok(result);
  }

  /**
   * \brief Parse fbm parameters
   */
  static MeshHexer::Result<FbmMeshParameters, std::string> parse_fbm_args(int& argc, char*** argv)
  {
    using Result = MeshHexer::Result<FbmMeshParameters, std::string>;

    FbmMeshParameters result;

    while(argc > 0)
    {
      char* cmd = (*argv)[0];

      if(cmp_argument("--help", cmd) || cmp_argument("-h", cmd))
      {
        consume_arg(argc, argv);
        result.show_help = true;
      }
      else if(cmp_argument("--output", cmd))
      {
        result.output = parse_argument("--output", argc, argv);
      }
      else if(cmp_argument("--level", cmd))
      {
        std::string level(parse_argument("--level", argc, argv));
        try
        {
          result.levels = std::stoull(level);
        }
        catch(const std::exception& e)
        {
          return Result::err(e.what());
        }
      }
      else if(cmp_argument("--bounding-box", cmd))
      {
        std::vector<char*> bb_args = parse_arguments("--bounding-box", argc, argv, 6);

        try
        {
          MeshHexer::BoundingBox bb{};
          bb.min.x = std::stod(bb_args[0]);
          bb.min.y = std::stod(bb_args[1]);
          bb.min.z = std::stod(bb_args[2]);
          bb.max.x = std::stod(bb_args[3]);
          bb.max.y = std::stod(bb_args[4]);
          bb.max.z = std::stod(bb_args[5]);

          result.bounding_box = bb;
        }
        catch(const std::exception& e)
        {
          return Result::err(e.what());
        }
      }
      else
      {
        // Unknown argument
        break;
      }
    }

    if(argc > 0)
    {
      result.mesh_file = (*argv[0]);
      consume_arg(argc, argv);
    }
    else if(!result.show_help)
    {
      return Result::err("Expected mesh file!");
    }

    return Result::ok(result);
  }

  /**
   * \brief Parse report parameters
   */
  static MeshHexer::Result<ReportParameters, std::string> parse_report_args(int& argc, char*** argv)
  {
    using Result = MeshHexer::Result<ReportParameters, std::string>;

    ReportParameters result;

    while(argc > 0)
    {
      char* cmd = (*argv)[0];

      if(cmp_argument("--help", cmd) || cmp_argument("-h", cmd))
      {
        consume_arg(argc, argv);
        result.show_help = true;
      }
      else if(cmp_argument("--configure-case-for-preprocessing", cmd))
      {
        consume_arg(argc, argv);
        result.configure_case_for_preprocessing = true;
      }
      else
      {
        // Unknown argument
        break;
      }
    }

    if(argc > 0)
    {
      result.mesh_file = (*argv[0]);
      consume_arg(argc, argv);
    }
    else if(!result.show_help)
    {
      return Result::err("Expected mesh file!");
    }

    return Result::ok(result);
  }

  /**
   * \brief Parse report parameters
   */
  static MeshHexer::Result<WarningsParameters, std::string> parse_warnings_args(int& argc, char*** argv)
  {
    using Result = MeshHexer::Result<WarningsParameters, std::string>;

    WarningsParameters result;

    while(argc > 0)
    {
      char* cmd = (*argv)[0];

      if(cmp_argument("--help", cmd) || cmp_argument("-h", cmd))
      {
        consume_arg(argc, argv);
        result.show_help = true;
      }
      if(cmp_argument("--summarize", cmd))
      {
        consume_arg(argc, argv);
        result.summarize = true;
      }
      else
      {
        // Unknown argument
        break;
      }
    }

    if(argc > 0)
    {
      result.mesh_file = (*argv[0]);
      consume_arg(argc, argv);
    }
    else if(!result.show_help)
    {
      return Result::err("Expected mesh file!");
    }

    return Result::ok(result);
  }

  static MeshHexer::Result<MeshHexer::Mesh, std::string> load_editable_mesh(const std::filesystem::path& filename)
  {
    using Result = MeshHexer::Result<MeshHexer::Mesh, std::string>;

    MeshHexer::Mesh mesh;
    std::string extension = filename.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) {
      return static_cast<char>(std::tolower(c));
    });

    if(extension == ".ply")
    {
      std::ifstream mesh_file(filename);
      if(!mesh_file)
      {
        return Result::err("Failed to open mesh " + filename.string());
      }
      std::string comment;
      if(!CGAL::IO::read_PLY(mesh_file, mesh, comment, true))
      {
        return Result::err("Failed to read mesh " + filename.string());
      }
    }
    else
    {
      if(!PMP::IO::read_polygon_mesh(filename.string(), mesh))
      {
        return Result::err("Failed to read mesh " + filename.string());
      }
    }

    if(CGAL::is_empty(mesh))
    {
      return Result::err("Mesh " + filename.string() + " is empty.");
    }

    if(!CGAL::is_triangle_mesh(mesh))
    {
      PMP::triangulate_faces(mesh);
    }

    return Result::ok(std::move(mesh));
  }

  int main(int argc, char* argv[])
  {
    int orig_argc = argc;
    char** orig_argv = argv;

    // Skip binary name
    consume_arg(argc, &argv);

    MeshHexer::Result<GlobalParameters, std::string> global_parse_result = parse_global_args(argc, &argv);

    if(global_parse_result.is_err())
    {
      std::cerr << "Parameter parsing failed with: " << global_parse_result.err_ref() << "\n";
      exit(1);
    }

    GlobalParameters gparams = global_parse_result.ok_value();

    if(gparams.show_help)
    {
      std::cout << usage;
      return 0;
    }

    if(gparams.command == "fbm-mesh")
    {
      MeshHexer::Result<FbmMeshParameters, std::string> parse_result = parse_fbm_args(argc, &argv);

      if(parse_result.is_err())
      {
        std::cerr << "Parameter parsing for command 'fbm-mesh' failed with: " << parse_result.err_ref() << "\n";
        exit(1);
      }

      FbmMeshParameters params = parse_result.ok_value();

      if(params.show_help)
      {
        std::cout << fbm_usage;
        exit(0);
      }

      MeshHexer::Result<MeshHexer::SurfaceMesh, std::string> result = MeshHexer::load_from_file(params.mesh_file, true);
      if(result.is_err())
      {
        std::cout << "Reading mesh failed with error: " << result.err_ref() << "\n";
        exit(1);
      }

      MeshHexer::SurfaceMesh mesh = std::move(result).take_ok();

      MeshHexer::FBMMeshSettings settings;

      settings.bounding_box = params.bounding_box.value_or(mesh.bounding_box());
      settings.levels = params.levels;

      MeshHexer::VolumeMesh vmesh = mesh.fbm_mesh(settings);

      std::ofstream mesh_file(params.output + ".xml");

      if(mesh_file.fail())
      {
        std::cerr << "Error opening " << params.output << ".xml for writing\n";
        return 1;
      }

      vmesh.write_feat_xml(mesh_file);

      std::ofstream mtx_file(params.output + ".mtx");

      if(mtx_file.fail())
      {
        std::cerr << "Error opening " << params.output << ".mtx for writing\n";
      }

      std::vector<std::uint64_t> sdls;
      sdls.reserve(vmesh.num_vertices());
      for(std::size_t i(0); i < vmesh.num_vertices(); i++)
      {
        sdls.push_back(vmesh.subdivision_level(i));
      }

      write_range_as_mtx(mtx_file, sdls.begin(), sdls.end());

      if(!gparams.checkpoint_path.empty())
      {
        MeshHexer::Result<void, std::string> result = mesh.write_to_file(gparams.checkpoint_path);

        if(result.is_err())
        {
          std::cout << "Writing checkpoint failed with error: " << result.err_ref() << "\n";
        }
      }
    }

    if(gparams.command == "min-gap")
    {
      MeshHexer::Result<MinGapParameters, std::string> parse_result = parse_mingap_args(argc, &argv);

      if(parse_result.is_err())
      {
        std::cerr << "Parameter parsing for command 'min-gap' failed with: " << parse_result.err_ref() << "\n";
        exit(1);
      }

      MinGapParameters params = parse_result.ok_value();

      if(params.show_help)
      {
        std::cout << mingap_usage;
        exit(0);
      }

      MeshHexer::Result<MeshHexer::SurfaceMesh, std::string> mesh_load_result = MeshHexer::load_from_file(params.mesh_file, true);
      if(mesh_load_result.is_err())
      {
        std::cout << "Reading mesh failed with error: " << mesh_load_result.err_ref() << "\n";
        exit(1);
      }

      MeshHexer::SurfaceMesh mesh = std::move(mesh_load_result).take_ok();

      const std::vector<std::size_t> short_normal_faces = mesh.faces_with_short_normals(1.0e-12);

      MeshHexer::Result<MeshHexer::Mesh, std::string> editable_mesh_result = load_editable_mesh(params.mesh_file);
      if(editable_mesh_result.is_err())
      {
        std::cerr << "Reading editable mesh failed with error: " << editable_mesh_result.err_ref() << "\n";
        exit(1);
      }

      MeshHexer::Mesh editable_mesh = std::move(editable_mesh_result).take_ok();
      const std::vector<MeshHexer::FaceIndex> small_area_face_indices =
        MeshHexer::faces_with_small_area(editable_mesh, SMALL_FACE_RELATIVE_AREA);

      std::unordered_set<std::size_t> excluded_face_ids;
      excluded_face_ids.reserve(short_normal_faces.size() + small_area_face_indices.size());
      for(std::size_t face_id : short_normal_faces)
      {
        excluded_face_ids.insert(face_id);
      }
      for(MeshHexer::FaceIndex face_index : small_area_face_indices)
      {
        excluded_face_ids.insert(face_index.idx());
      }

      std::cerr << "Min-gap analysis filter: "
                << small_area_face_indices.size() << " near-zero-area faces, "
                << short_normal_faces.size() << " faces with collapsed centroid normals";
      if(excluded_face_ids.empty())
      {
        std::cerr << " excluded from analysis: 0\n";
      }
      else
      {
        std::cerr << " excluded from analysis: " << excluded_face_ids.size() << "\n";
      }

      std::filesystem::path temp_analysis_mesh_path;
      if(!excluded_face_ids.empty())
      {
        for(std::size_t face_id : excluded_face_ids)
        {
          editable_mesh.remove_face(MeshHexer::FaceIndex(face_id));
        }
        editable_mesh.collect_garbage();

        temp_analysis_mesh_path =
          std::filesystem::temp_directory_path() /
          (params.mesh_file.stem().string() + "_mingap_filtered_" + std::to_string(getpid()) + ".off");

        if(!CGAL::IO::write_polygon_mesh(temp_analysis_mesh_path.string(), editable_mesh))
        {
          std::cerr << "Failed to write filtered analysis mesh to " << temp_analysis_mesh_path << "\n";
          exit(1);
        }

        MeshHexer::Result<MeshHexer::SurfaceMesh, std::string> filtered_mesh_result =
          MeshHexer::load_from_file(temp_analysis_mesh_path, true);
        if(filtered_mesh_result.is_err())
        {
          std::filesystem::remove(temp_analysis_mesh_path);
          std::cerr << "Reading filtered analysis mesh failed with error: " << filtered_mesh_result.err_ref() << "\n";
          exit(1);
        }
        mesh = std::move(filtered_mesh_result).take_ok();
      }

      MeshHexer::Gap min_gap = mesh.min_gap();
      const double coarse_mesh_size = min_gap.diameter * MeshHexer::min_gap_config().coarse_mesh_scaling *
        std::pow(3.0, static_cast<double>(min_gap.bin_span));
      std::cerr << "Min-gap monitor histogram (area-based, factor-of-3 bins; Collections = connected agglomerates "
                   "with area >= "
                << (MeshHexer::min_gap_config().histogram_min_fraction_connected_threshold * 100.0)
                << "% of total surface using faces from this bin only; Connected/Isolated/Max Component % "
                   "are measured against the total histogram-covered area)\n";
      mesh.write_monitor_histogram(std::cerr);
      print_min_gap(min_gap, params.verbose, coarse_mesh_size);

      if(!temp_analysis_mesh_path.empty())
      {
        std::filesystem::remove(temp_analysis_mesh_path);
      }

      if(!gparams.checkpoint_path.empty())
      {
        MeshHexer::Result<void, std::string> result = mesh.write_to_file(gparams.checkpoint_path);

        if(result.is_err())
        {
          std::cout << "Writing checkpoint failed with error: " << result.err_ref() << "\n";
        }
      }
    }

    if(gparams.command == "report")
    {
      MeshHexer::Result<ReportParameters, std::string> parse_result = parse_report_args(argc, &argv);

      if(parse_result.is_err())
      {
        std::cerr << "Parameter parsing for command 'report' failed with: " << parse_result.err_ref() << "\n";
        exit(1);
      }

      ReportParameters params = parse_result.ok_value();

      if(params.show_help)
      {
        std::cout << report_usage;
        exit(0);
      }

      if(params.configure_case_for_preprocessing)
      {
        MeshHexer::Result<std::filesystem::path, std::string> case_dir_result =
          resolve_case_directory(params.mesh_file);
        if(case_dir_result.is_err())
        {
          std::cerr << "Case input resolution failed with: " << case_dir_result.err_ref() << "\n";
          exit(1);
        }

        const std::filesystem::path case_dir = std::filesystem::canonical(case_dir_result.ok_value());
        const std::filesystem::path setup_path = case_dir / "setup.e3d";
        MeshHexer::Result<ParsedRoundSetup, std::string> setup_result = parse_round_setup_file(setup_path);
        if(setup_result.is_err())
        {
          std::cerr << "Round geometry setup parsing failed with: " << setup_result.err_ref() << "\n";
          exit(1);
        }

        const MeshHexer::RoundGeometryAnalysisConfig& analysis_config = MeshHexer::round_geometry_analysis_config();
        ParsedRoundSetup parsed_setup = setup_result.ok_value();
        std::string geometry_type = parsed_setup.geometry_type;
        bool geometry_type_from_default = false;
        if(!parsed_setup.has_geometry_type)
        {
          geometry_type = analysis_config.default_geometry_type;
          geometry_type_from_default = true;
        }

        if(geometry_type == "BOX")
        {
          MeshHexer::Result<std::filesystem::path, std::string> mesh_path_result =
            resolve_report_mesh_path(case_dir);
          if(mesh_path_result.is_err())
          {
            std::cerr << "Report input resolution failed with: " << mesh_path_result.err_ref() << "\n";
            exit(1);
          }

          const std::filesystem::path mesh_path = mesh_path_result.ok_value();
          MeshHexer::Result<MeshHexer::SurfaceMesh, std::string> result = MeshHexer::load_from_file(mesh_path, true);
          if(result.is_err())
          {
            std::cout << "Reading mesh failed with error: " << result.err_ref() << "\n";
            exit(1);
          }

          MeshHexer::SurfaceMesh mesh = std::move(result).take_ok();
          std::filesystem::path absolute_path = std::filesystem::canonical(mesh_path);
          MeshHexer::MeshWarnings warnings = mesh.warnings();

          const bool box_needs_population = parsed_setup.has_geometry_start && parsed_setup.has_geometry_length;
          const bool box_already_prepared = has_preconfigured_box_preprocessing(parsed_setup);
          bool setup_rewritten = false;
          std::filesystem::path backup_setup_path = case_dir / "setup_BU.e3d";

          if(box_needs_population)
          {
            MeshHexer::Result<SetupRewriteResult, std::string> rewrite_result =
              rewrite_setup_for_box_geometry(setup_path, parsed_setup, analysis_config);
            if(rewrite_result.is_err())
            {
              std::cerr << "BOX setup rewrite failed with: " << rewrite_result.err_ref() << "\n";
              exit(1);
            }

            const std::filesystem::path temp_setup_path = case_dir / "setup.e3d.tmp";
            std::error_code ec;
            std::filesystem::copy_file(
              setup_path,
              backup_setup_path,
              std::filesystem::copy_options::overwrite_existing,
              ec);
            if(ec)
            {
              std::cerr << "Failed to back up setup file to " << backup_setup_path << ": " << ec.message() << "\n";
              exit(1);
            }
            if(!write_text_lines(temp_setup_path, rewrite_result.ok_ref().lines))
            {
              std::cerr << "Failed to write updated setup file " << temp_setup_path << "\n";
              std::filesystem::remove(temp_setup_path, ec);
              exit(1);
            }
            std::filesystem::rename(temp_setup_path, setup_path, ec);
            if(ec)
            {
              std::cerr << "Failed to replace setup file " << setup_path << ": " << ec.message() << "\n";
              std::filesystem::remove(temp_setup_path, ec);
              exit(1);
            }
            setup_rewritten = true;
          }
          else if(!box_already_prepared)
          {
            std::cerr << "BOX setup rewrite failed: neither geometryStart/geometryLength in [E3DGeometryData/Machine] "
                      << "nor a prepared [E3DGeometryData/Preprocessing] Box section were found in "
                      << setup_path << "\n";
            exit(1);
          }

          std::cout << Markdown::h1("Mesh-Report for " + absolute_path.filename().string()) << "\n\n";
          std::cout << Markdown::h2("Metadata") << "\n";
          std::cout << Markdown::li("Path: " + absolute_path.string()) << "\n\n";
          std::cout << Markdown::h2("Topology") << "\n";
          std::cout << Markdown::li("Number of vertices: " + std::to_string(mesh.num_vertices())) << "\n";
          std::cout << Markdown::li("Number of edges: " + std::to_string(mesh.num_edges())) << "\n";
          std::cout << Markdown::li("Number of faces: " + std::to_string(mesh.num_faces())) << "\n";
          MeshHexer::BoundingBox bb = mesh.bounding_box();
          std::cout << Markdown::li("Extent (x, y, z): [" +
            std::to_string(bb.min.x) + ", " +
            std::to_string(bb.max.x) + "] x [" +
            std::to_string(bb.min.y) + ", " +
            std::to_string(bb.max.y) + "] x [" +
            std::to_string(bb.min.z) + ", " +
            std::to_string(bb.max.z) + "]") << "\n";
          std::cout << Markdown::li("Is closed: " + std::string(mesh.is_closed() ? "True" : "False")) << "\n";
          std::cout << Markdown::li(
                         "Is wound consistently: " + std::string(mesh.is_wound_consistently() ? "True" : "False"))
                    << "\n";
          std::cout << Markdown::li("Is oriented outward: " + std::string(mesh.is_outward_oriented() ? "True" : "False"))
                    << "\n";
          std::cout << Markdown::li("Minimal triangle aspect ratio: " + std::to_string(mesh.minimal_aspect_ratio()))
                    << "\n";
          std::cout << Markdown::li("Maximal triangle aspect ratio: " + std::to_string(mesh.maximal_aspect_ratio()))
                    << "\n\n";
          std::cout << Markdown::h2("Defects") << "\n";
          std::cout << Markdown::li("Self-intersections: " + std::to_string(warnings.self_intersections.size())) << "\n";
          std::cout << Markdown::li("Degenerate triangles: " + std::to_string(warnings.degenerate_triangles.size()))
                    << "\n";
          std::cout << Markdown::li("Anisotropic triangles: " + std::to_string(warnings.anisotropic_triangles.size()))
                    << "\n\n";
          std::cout << Markdown::h2("Box Preprocessing") << "\n";
          std::cout << Markdown::li("setup.e3d: " + setup_path.string()) << "\n";
          if(geometry_type_from_default)
          {
            std::cout << Markdown::li("GeometryType: defaulted to " + geometry_type + " from preprocessor.cfg") << "\n";
          }
          if(setup_rewritten)
          {
            std::cout << Markdown::li("setup backup: " + backup_setup_path.string()) << "\n";
            std::cout << Markdown::li("geometryStart moved to preprocessing: " + parsed_setup.geometry_start_raw) << "\n";
            std::cout << Markdown::li("geometryLength moved to preprocessing: " + parsed_setup.geometry_length_raw) << "\n";
            std::cout << Markdown::li("setup.e3d preprocessing section refreshed for Box") << "\n\n";
          }
          else
          {
            std::cout << Markdown::li(
                           "setup.e3d preprocessing section already prepared for Box; setup file left unchanged")
                      << "\n\n";
          }
          return 0;
        }
        if(geometry_type != "ROUND")
        {
          std::cerr << "Case configuration failed: unsupported GeometryType=" << geometry_type
                    << "\n";
          exit(1);
        }

        MeshHexer::Result<std::filesystem::path, std::string> mesh_path_result =
          resolve_report_mesh_path(case_dir);
        if(mesh_path_result.is_err())
        {
          std::cerr << "Report input resolution failed with: " << mesh_path_result.err_ref() << "\n";
          exit(1);
        }

        const std::filesystem::path mesh_path = mesh_path_result.ok_value();
        MeshHexer::Result<MeshHexer::SurfaceMesh, std::string> result = MeshHexer::load_from_file(mesh_path, true);
        if(result.is_err())
        {
          std::cout << "Reading mesh failed with error: " << result.err_ref() << "\n";
          exit(1);
        }

        MeshHexer::SurfaceMesh mesh = std::move(result).take_ok();
        std::filesystem::path absolute_path = std::filesystem::canonical(mesh_path);
        MeshHexer::MeshWarnings warnings = mesh.warnings();

        std::cout << Markdown::h1("Mesh-Report for " + absolute_path.filename().string()) << "\n\n";
        std::cout << Markdown::h2("Metadata") << "\n";
        std::cout << Markdown::li("Path: " + absolute_path.string()) << "\n\n";
        std::cout << Markdown::h2("Topology") << "\n";
        std::cout << Markdown::li("Number of vertices: " + std::to_string(mesh.num_vertices())) << "\n";
        std::cout << Markdown::li("Number of edges: " + std::to_string(mesh.num_edges())) << "\n";
        std::cout << Markdown::li("Number of faces: " + std::to_string(mesh.num_faces())) << "\n";
        MeshHexer::BoundingBox bb = mesh.bounding_box();
        std::cout << Markdown::li("Extent (x, y, z): [" +
          std::to_string(bb.min.x) + ", " +
          std::to_string(bb.max.x) + "] x [" +
          std::to_string(bb.min.y) + ", " +
          std::to_string(bb.max.y) + "] x [" +
          std::to_string(bb.min.z) + ", " +
          std::to_string(bb.max.z) + "]") << "\n";
        std::cout << Markdown::li("Is closed: " + std::string(mesh.is_closed() ? "True" : "False")) << "\n";
        std::cout << Markdown::li(
                       "Is wound consistently: " + std::string(mesh.is_wound_consistently() ? "True" : "False"))
                  << "\n";
        std::cout << Markdown::li("Is oriented outward: " + std::string(mesh.is_outward_oriented() ? "True" : "False"))
                  << "\n";
        std::cout << Markdown::li("Minimal triangle aspect ratio: " + std::to_string(mesh.minimal_aspect_ratio()))
                  << "\n";
        std::cout << Markdown::li("Maximal triangle aspect ratio: " + std::to_string(mesh.maximal_aspect_ratio()))
                  << "\n\n";
        std::cout << Markdown::h2("Defects") << "\n";
        std::cout << Markdown::li("Self-intersections: " + std::to_string(warnings.self_intersections.size())) << "\n";
        std::cout << Markdown::li("Degenerate triangles: " + std::to_string(warnings.degenerate_triangles.size()))
                  << "\n";
        std::cout << Markdown::li("Anisotropic triangles: " + std::to_string(warnings.anisotropic_triangles.size()))
                  << "\n\n";

        MeshHexer::Result<MeshHexer::RoundGeometryAnalysisResult, std::string> analysis_result =
          mesh.round_geometry_analysis(parsed_setup.inflows, analysis_config);
        if(analysis_result.is_err())
        {
          std::cerr << "Round geometry analysis failed with: " << analysis_result.err_ref() << "\n";
          exit(1);
        }

        const MeshHexer::RoundGeometryAnalysisResult& analysis = analysis_result.ok_ref();

        std::cout << Markdown::h2("Round Geometry Analysis") << "\n";
        std::cout << Markdown::li("setup.e3d: " + setup_path.string()) << "\n";
        std::cout << Markdown::li(
                       "Axis center (x, y): (" + std::to_string(analysis.axis_center.x) + ", " +
                       std::to_string(analysis.axis_center.y) + ")")
                  << "\n";
        std::cout << Markdown::li(
                       "Axis aligned to origin: " + std::string(analysis.axis_aligned_to_origin ? "True" : "False"))
                  << "\n";
        std::cout << Markdown::li("Top slice z: " + std::to_string(analysis.top_slice_z)) << "\n";
        std::cout << Markdown::li("Outer diameter: " + std::to_string(analysis.outer_diameter)) << "\n";
        std::cout << Markdown::li("Inner diameter: " + std::to_string(analysis.inner_diameter)) << "\n";
        std::cout << Markdown::li("Inner/outer ratio: " + std::to_string(analysis.inner_to_outer_ratio)) << "\n";
        std::cout << Markdown::li("Classification: " + classification_to_string(analysis.classification)) << "\n";
        std::cout << Markdown::li("Extrusion length: " + std::to_string(analysis.extrusion_length)) << "\n";
        std::cout << Markdown::li(
                       "Extrusion length consistent: " +
                       std::string(analysis.extrusion_length_consistent ? "True" : "False"))
                  << "\n";
        std::cout << Markdown::li(
                       "z-min limited by axial inflow: " +
                       std::string(analysis.z_min_limited_by_axial_inflow ? "True" : "False"))
                  << "\n";
        std::cout << Markdown::li("z_min_phys: " + std::to_string(analysis.z_min_physical)) << "\n";
        std::cout << Markdown::li("z_max_phys: " + std::to_string(analysis.z_max_physical)) << "\n";
        if(!analysis.extrusion_length_samples.empty())
        {
          std::ostringstream sample_stream;
          for(std::size_t i = 0; i < analysis.extrusion_length_samples.size(); ++i)
          {
            if(i > 0)
            {
              sample_stream << ", ";
            }
            sample_stream << analysis.extrusion_length_samples[i];
          }
          std::cout << Markdown::li("Extrusion samples: [" + sample_stream.str() + "]") << "\n";
        }

        if(params.configure_case_for_preprocessing)
        {
          const double shift_x = analysis.axis_aligned_to_origin ? 0.0 : -analysis.axis_center.x;
          const double shift_y = analysis.axis_aligned_to_origin ? 0.0 : -analysis.axis_center.y;
          const bool setup_preprocessing_matches =
            round_preprocessing_section_matches(setup_path, analysis, analysis_config);
          const bool setup_requires_update = !analysis.axis_aligned_to_origin || !setup_preprocessing_matches;

          std::optional<SetupRewriteResult> rewrite_result;
          const std::filesystem::path backup_surface_path = absolute_path.parent_path() / "surface_BU.off";
          const std::filesystem::path backup_setup_path = absolute_path.parent_path() / "setup_BU.e3d";
          const std::string mesh_extension = absolute_path.has_extension() ? absolute_path.extension().string() : ".off";
          const std::filesystem::path temp_surface_path =
            absolute_path.parent_path() / (absolute_path.stem().string() + ".tmp" + mesh_extension);
          const std::filesystem::path temp_setup_path = absolute_path.parent_path() / "setup.e3d.tmp";

          if(setup_requires_update)
          {
            MeshHexer::Result<SetupRewriteResult, std::string> rewrite_attempt =
              rewrite_setup_for_round_geometry(setup_path, shift_x, shift_y, analysis, analysis_config);
            if(rewrite_attempt.is_err())
            {
              std::cerr << "Round geometry setup rewrite failed with: " << rewrite_attempt.err_ref() << "\n";
              exit(1);
            }
            rewrite_result = rewrite_attempt.ok_value();

            std::error_code ec;
            std::filesystem::copy_file(
              setup_path,
              backup_setup_path,
              std::filesystem::copy_options::overwrite_existing,
              ec);
            if(ec)
            {
              std::cerr << "Failed to back up setup file to " << backup_setup_path << ": " << ec.message() << "\n";
              exit(1);
            }

            if(!analysis.axis_aligned_to_origin)
            {
              std::filesystem::copy_file(
                absolute_path,
                backup_surface_path,
                std::filesystem::copy_options::overwrite_existing,
                ec);
              if(ec)
              {
                std::cerr << "Failed to back up surface mesh to " << backup_surface_path << ": " << ec.message() << "\n";
                exit(1);
              }

              MeshHexer::SurfaceMesh shifted_mesh = std::move(mesh);
              shifted_mesh.translate(shift_x, shift_y, 0.0);

              MeshHexer::Result<void, std::string> write_result = shifted_mesh.write_to_file(temp_surface_path.string());
              if(write_result.is_err())
              {
                std::cerr << "Failed to write shifted mesh: " << write_result.err_ref() << "\n";
                std::filesystem::remove(temp_surface_path, ec);
                exit(1);
              }
            }

            if(!write_text_lines(temp_setup_path, rewrite_result->lines))
            {
              std::cerr << "Failed to write updated setup file " << temp_setup_path << "\n";
              std::filesystem::remove(temp_surface_path, ec);
              std::filesystem::remove(temp_setup_path, ec);
              exit(1);
            }

            if(!analysis.axis_aligned_to_origin)
            {
              std::filesystem::rename(temp_surface_path, absolute_path, ec);
              if(ec)
              {
                std::cerr << "Failed to replace surface mesh " << absolute_path << ": " << ec.message() << "\n";
                std::filesystem::remove(temp_surface_path, ec);
                std::filesystem::remove(temp_setup_path, ec);
                exit(1);
              }
            }

            std::filesystem::rename(temp_setup_path, setup_path, ec);
            if(ec)
            {
              std::cerr << "Failed to replace setup file " << setup_path << ": " << ec.message() << "\n";
              std::filesystem::remove(temp_setup_path, ec);
              exit(1);
            }
          }

          if(analysis.axis_aligned_to_origin)
          {
            std::cout << Markdown::li("Axis fix: not required") << "\n";
          }
          else
          {
            std::cout << Markdown::li(
                           "Axis fix: applied shift (" + std::to_string(shift_x) + ", " + std::to_string(shift_y) +
                           ", 0.000000)")
                      << "\n";
            std::cout << Markdown::li("surface backup: " + backup_surface_path.string()) << "\n";
          }
          if(!setup_requires_update)
          {
            std::cout << Markdown::li("setup.e3d preprocessing section already correct; setup file left unchanged") << "\n";
          }
          else
          {
            std::cout << Markdown::li("setup backup: " + backup_setup_path.string()) << "\n";
            std::cout << Markdown::li(
                           "setup.e3d updates: shifted " + std::to_string(rewrite_result->shifted_centers) +
                           " center entries and " + std::to_string(rewrite_result->shifted_midpoints) +
                           " midpoint entries; preprocessing section refreshed")
                      << "\n";
          }
        }

        std::cout << "\n";
      }
    }

    if(gparams.command == "warnings")
    {
      MeshHexer::Result<WarningsParameters, std::string> parse_result = parse_warnings_args(argc, &argv);

      if(parse_result.is_err())
      {
        std::cerr << "Parameter parsing for command 'warnings' failed with: " << parse_result.err_ref() << "\n";
        exit(1);
      }

      WarningsParameters params = parse_result.ok_value();

      if(params.show_help)
      {
        std::cout << warnings_usage;
        exit(0);
      }

      MeshHexer::Result<MeshHexer::SurfaceMesh, std::string> result = MeshHexer::load_from_file(params.mesh_file, true);
      if(result.is_err())
      {
        std::cout << "Reading mesh failed with error: " << result.err_ref() << "\n";
        exit(1);
      }

      MeshHexer::SurfaceMesh mesh = std::move(result).take_ok();
      MeshHexer::MeshWarnings warnings = mesh.warnings();

      if(params.summarize)
      {
        std::cout << std::to_string(warnings.self_intersections.size()) << " x Self-intersection of mesh ["
                  << MeshHexer::SelfIntersectionWarning::name << "]\n";
      }
      else
      {
        for(MeshHexer::SelfIntersectionWarning& warning : warnings.self_intersections)
        {
          std::cout << "Self-intersection of mesh between triangle " << std::to_string(warning.tri_a)
                    << " and triangle " << std::to_string(warning.tri_b) << " ["
                    << MeshHexer::SelfIntersectionWarning::name << "]\n";
        }
      }

      if(params.summarize)
      {
        std::cout << std::to_string(warnings.degenerate_triangles.size()) << " x Triangle with colinear coordinates ["
                  << MeshHexer::DegenerateTriangleWarning::name << "]\n";
      }
      else
      {
        for(MeshHexer::DegenerateTriangleWarning& warning : warnings.degenerate_triangles)
        {
          std::cout << "Coordinates of triangle " << std::to_string(warning.idx) << " are colinear ["
                    << MeshHexer::DegenerateTriangleWarning::name << "]\n";
        }
      }

      if(params.summarize)
      {
        std::cout << std::to_string(warnings.anisotropic_triangles.size())
                  << " x Triangle with aspect ratio greater than 15 [" << MeshHexer::DegenerateTriangleWarning::name
                  << "]\n";
      }
      else
      {
        for(MeshHexer::AnisotropicTriangleWarning& warning : warnings.anisotropic_triangles)
        {
          std::cout << "Triangle " << std::to_string(warning.idx) << " has aspect ratio greater than 15 ["
                    << MeshHexer::DegenerateTriangleWarning::name << "]\n";
        }
      }
    }

    return 0;
  }
} // namespace MeshHexerCLI

int main(int argc, char* argv[])
{
  MeshHexerCLI::main(argc, argv);
}
