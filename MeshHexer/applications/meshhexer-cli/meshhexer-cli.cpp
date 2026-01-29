#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <meshhexer/meshhexer.hpp>
#include <meshhexer/types.hpp>

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
    void print_min_gap(const MeshHexer::Gap& min_gap, bool verbose)
    {
      if(verbose)
      {
        std::cout << "Min-gap of " << min_gap.diameter << " between faces " << min_gap.face << " and " << min_gap.opposite_face << "\n";
        std::cout << "Use `SelectIDs(IDs=[0, " << min_gap.face << ", 0, " << min_gap.opposite_face
              << "], FieldType='CELL')` to select the chosen triangles in ParaView\n";
      }
      else
      {
        std::cout << min_gap.diameter << "\n";
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
    "Usage: meshhexer-cli report <mesh>\n"
    "\n"
    "Print information about the mesh\n";

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

      MeshHexer::Result<MeshHexer::SurfaceMesh, std::string> result = MeshHexer::load_from_file(params.mesh_file, true);
      if(result.is_err())
      {
        std::cout << "Reading mesh failed with error: " << result.err_ref() << "\n";
        exit(1);
      }

      MeshHexer::SurfaceMesh mesh = std::move(result).take_ok();

      MeshHexer::Gap min_gap = mesh.min_gap();
      print_min_gap(min_gap, params.verbose);

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

      MeshHexer::Result<MeshHexer::SurfaceMesh, std::string> result = MeshHexer::load_from_file(params.mesh_file, true);
      if(result.is_err())
      {
        std::cout << "Reading mesh failed with error: " << result.err_ref() << "\n";
        exit(1);
      }

      MeshHexer::SurfaceMesh mesh = std::move(result).take_ok();

      std::filesystem::path absolute_path = std::filesystem::canonical(params.mesh_file);

      MeshHexer::MeshWarnings warnings = mesh.warnings();

      std::cout << Markdown::h1("Mesh-Report for " + absolute_path.filename().string()) << "\n";

      std::cout << "\n";

      std::cout << Markdown::h2("Metadata") << "\n";
      std::cout << Markdown::li("Path: " + absolute_path.string()) << "\n";

      std::cout << "\n";

      std::cout << Markdown::h2("Topology") << "\n";
      std::cout << Markdown::li("Number of vertices: " + std::to_string(mesh.num_vertices())) << "\n";
      std::cout << Markdown::li("Number of edges: " + std::to_string(mesh.num_edges())) << "\n";
      std::cout << Markdown::li("Number of faces: " + std::to_string(mesh.num_faces())) << "\n";

      MeshHexer::BoundingBox bb = mesh.bounding_box();
      // clang-format off
      std::cout << Markdown::li("Extent (x, y, z): [" +
        std::to_string(bb.min.x) + ", " +
        std::to_string(bb.max.x) + "] x [" +
        std::to_string(bb.min.y) + ", " +
        std::to_string(bb.max.y) + "] x [" +
        std::to_string(bb.min.z) + ", " +
        std::to_string(bb.max.z) + "]") << "\n";
      // clang-format on

      std::cout << Markdown::li("Is closed: " + std::string(mesh.is_closed() ? "True" : "False")) << "\n";
      std::cout << Markdown::li(
                     "Is wound consistently: " + std::string(mesh.is_wound_consistently() ? "True" : "False"))
                << "\n";
      std::cout << Markdown::li("Is oriented outward: " + std::string(mesh.is_outward_oriented() ? "True" : "False"))
                << "\n";
      std::cout << Markdown::li("Minimal triangle aspect ratio: " + std::to_string(mesh.minimal_aspect_ratio()))
                << "\n";
      std::cout << Markdown::li("Maximal triangle aspect ratio: " + std::to_string(mesh.maximal_aspect_ratio()))
                << "\n";

      std::cout << "\n";

      std::cout << Markdown::h2("Defects") << "\n";
      std::cout << Markdown::li("Self-intersections: " + std::to_string(warnings.self_intersections.size())) << "\n";
      std::cout << Markdown::li("Degenerate triangles: " + std::to_string(warnings.degenerate_triangles.size()))
                << "\n";
      std::cout << Markdown::li("Anisotropic triangles: " + std::to_string(warnings.anisotropic_triangles.size()))
                << "\n";

      std::cout << "\n";
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
