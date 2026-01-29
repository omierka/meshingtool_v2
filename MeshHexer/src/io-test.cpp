#include <io.hpp>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <catch2/catch_test_macros.hpp>
#include <cgal_types.hpp>
#include <meshhexer/meshhexer.hpp>
#include <meshhexer/types.hpp>

#include <sstream>
#include <string_view>

using namespace MeshHexer;

namespace
{
  /**
   * \brief Returns a std::stringstream containing a OFF file of a unit cube
   */
  std::stringstream cube_mesh()
  {
    std::stringstream mts;
    mts << "OFF\n";
    mts << "8 12 0\n";
    mts << "-0.5 -0.5 -0.5\n";
    mts << "0.5 -0.5 -0.5\n";
    mts << "-0.5 0.5 -0.5\n";
    mts << "0.5 0.5 -0.5\n";
    mts << "-0.5 -0.5 0.5\n";
    mts << "0.5 -0.5 0.5\n";
    mts << "-0.5 0.5 0.5\n";
    mts << "0.5 0.5 0.5\n";
    mts << "3 1 0 3\n";
    mts << "3 0 2 3\n";
    mts << "3 4 5 6\n";
    mts << "3 5 7 6\n";
    mts << "3 0 1 4\n";
    mts << "3 1 5 4\n";
    mts << "3 3 2 7\n";
    mts << "3 2 6 7\n";
    mts << "3 2 0 6\n";
    mts << "3 0 4 6\n";
    mts << "3 1 3 5\n";
    mts << "3 3 7 5\n";

    return mts;
  }

  std::stringstream cube_vtu()
  {
    std::stringstream result;

    result << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    result << "<UnstructuredGrid>\n";
    result << "<Piece NumberOfPoints=\"8\" NumberOfCells=\"12\">\n";
    result << "<Points>\n";
    result << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    result << "-0.5 -0.5 -0.5\n";
    result << "0.5 -0.5 -0.5\n";
    result << "-0.5 0.5 -0.5\n";
    result << "0.5 0.5 -0.5\n";
    result << "-0.5 -0.5 0.5\n";
    result << "0.5 -0.5 0.5\n";
    result << "-0.5 0.5 0.5\n";
    result << "0.5 0.5 0.5\n";
    result << "</DataArray>\n";
    result << "</Points>\n";
    result << "<Cells>\n";
    result << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    result << "1\n";
    result << "0\n";
    result << "3\n";
    result << "0\n";
    result << "2\n";
    result << "3\n";
    result << "4\n";
    result << "5\n";
    result << "6\n";
    result << "5\n";
    result << "7\n";
    result << "6\n";
    result << "0\n";
    result << "1\n";
    result << "4\n";
    result << "1\n";
    result << "5\n";
    result << "4\n";
    result << "3\n";
    result << "2\n";
    result << "7\n";
    result << "2\n";
    result << "6\n";
    result << "7\n";
    result << "2\n";
    result << "0\n";
    result << "6\n";
    result << "0\n";
    result << "4\n";
    result << "6\n";
    result << "1\n";
    result << "3\n";
    result << "5\n";
    result << "3\n";
    result << "7\n";
    result << "5\n";
    result << "</DataArray>\n";
    result << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    result << "3\n";
    result << "6\n";
    result << "9\n";
    result << "12\n";
    result << "15\n";
    result << "18\n";
    result << "21\n";
    result << "24\n";
    result << "27\n";
    result << "30\n";
    result << "33\n";
    result << "36\n";
    result << "</DataArray>\n";
    result << "<DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    result << "5\n";
    result << "5\n";
    result << "5\n";
    result << "5\n";
    result << "5\n";
    result << "5\n";
    result << "5\n";
    result << "5\n";
    result << "5\n";
    result << "5\n";
    result << "5\n";
    result << "5\n";
    result << "</DataArray>\n";
    result << "</Cells>\n";
    result << "<CellData>\n";
    result << "<DataArray type=\"UInt32\" Name=\"data\" format=\"ascii\">\n";
    result << "0\n";
    result << "1\n";
    result << "2\n";
    result << "3\n";
    result << "4\n";
    result << "5\n";
    result << "6\n";
    result << "7\n";
    result << "8\n";
    result << "9\n";
    result << "10\n";
    result << "11\n";
    result << "</DataArray>\n";
    result << "</CellData>\n";
    result << "</Piece>\n";
    result << "</UnstructuredGrid>\n";
    result << "</VTKFile>\n";

    return result;
  }
} // namespace

TEST_CASE("Write FEAT mesh format", "[IO]")
{
  VolumeMesh vmesh(2, 2, 2);

  vmesh.xs_begin()[0] = 0.0;
  vmesh.xs_begin()[1] = 1.0;

  vmesh.ys_begin()[0] = 0.0;
  vmesh.ys_begin()[1] = 1.0;

  vmesh.zs_begin()[0] = 0.0;
  vmesh.zs_begin()[1] = 1.0;

  std::ostringstream output;

  vmesh.write_feat_xml(output);

  constexpr std::string_view expected("<FeatMeshFile version=\"1\" mesh=\"conformal:hypercube:3:3\">\n"
                                      "<Mesh type=\"conformal:hypercube:3:3\" size=\"8 12 6 1\">\n"
                                      "<Vertices>\n"
                                      "0 0 0\n"
                                      "1 0 0\n"
                                      "0 1 0\n"
                                      "1 1 0\n"
                                      "0 0 1\n"
                                      "1 0 1\n"
                                      "0 1 1\n"
                                      "1 1 1\n"
                                      "</Vertices>\n"
                                      "<Topology dim=\"1\">\n"
                                      "0 1\n"
                                      "2 3\n"
                                      "4 5\n"
                                      "6 7\n"
                                      "0 2\n"
                                      "1 3\n"
                                      "4 6\n"
                                      "5 7\n"
                                      "0 4\n"
                                      "1 5\n"
                                      "2 6\n"
                                      "3 7\n"
                                      "</Topology>\n"
                                      "<Topology dim=\"2\">\n"
                                      "0 1 2 3\n"
                                      "4 5 6 7\n"
                                      "0 1 4 5\n"
                                      "2 3 6 7\n"
                                      "0 2 4 6\n"
                                      "1 3 5 7\n"
                                      "</Topology>\n"
                                      "<Topology dim=\"3\">\n"
                                      "0 1 2 3 4 5 6 7\n"
                                      "</Topology>\n"
                                      "</Mesh>\n"
                                      "</FeatMeshFile>\n");

  REQUIRE(output.str() == expected);
}

TEST_CASE("Write mesh to .vtu", "[IO]")
{
  std::stringstream cube_file = cube_mesh();

  Mesh mesh;
  CGAL::IO::read_OFF(cube_file, mesh);

  auto prop = mesh.add_property_map<FaceIndex, std::uint32_t>("f:data");

  for(FaceIndex v : mesh.faces())
  {
    prop.first[v] = static_cast<std::uint32_t>(v);
  }

  std::stringstream vtu;
  write_vtu(vtu, mesh);

  REQUIRE(vtu.str() == cube_vtu().str());
}

TEST_CASE("Read mesh from VTU", "[IO]")
{
  std::stringstream cube_file = cube_vtu();
  Result<Mesh, std::string> result = read_vtu(cube_file);

  if(result.is_err())
  {
    std::cerr << result.err_ref() << "\n";
  }

  REQUIRE(!result.is_err());
  Mesh mesh = std::move(result).take_ok();

  REQUIRE(mesh.number_of_vertices() == 8);
  REQUIRE(mesh.number_of_edges() == 18);
  REQUIRE(mesh.number_of_faces() == 12);

  REQUIRE(mesh.point(VertexIndex(0)) == Point3D(-0.5, -0.5, -0.5));
  REQUIRE(mesh.point(VertexIndex(1)) == Point3D(0.5, -0.5, -0.5));
  REQUIRE(mesh.point(VertexIndex(2)) == Point3D(-0.5, 0.5, -0.5));
  REQUIRE(mesh.point(VertexIndex(3)) == Point3D(0.5, 0.5, -0.5));
  REQUIRE(mesh.point(VertexIndex(4)) == Point3D(-0.5, -0.5, 0.5));
  REQUIRE(mesh.point(VertexIndex(5)) == Point3D(0.5, -0.5, 0.5));
  REQUIRE(mesh.point(VertexIndex(6)) == Point3D(-0.5, 0.5, 0.5));
  REQUIRE(mesh.point(VertexIndex(7)) == Point3D(0.5, 0.5, 0.5));

  auto opt_prop_map = mesh.property_map<FaceIndex, std::uint32_t>("f:data");

  REQUIRE(opt_prop_map.has_value());

  for(std::size_t i(0); i < mesh.number_of_faces(); i++)
  {
    REQUIRE(opt_prop_map.value()[FaceIndex(i)] == i);
  }
}

TEST_CASE("Test XML lexing", "[IO]")
{
  std::stringstream xml;
  xml << "<?xml version=\"1.0\"?>\n";
  xml << "<!-- This is - a - test -->\n";
  xml << "<Foo test=\"1\">\n";
  xml << "<Bar test=\"2\">\n";
  xml << "String\n";
  xml << "<Baz/>\n";
  xml << "<Baz attr = \"Hello World!\" />\n";
  xml << "</Bar>\n";
  xml << "</Foo>\n";

  XML::XMLTokenizer tokenizer(xml);

  XML::Token token(XML::TokenType::TagBegin, "");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::DeclarationBegin);

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::AttributeName);
  REQUIRE(token.value == "version");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::Equals);

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::AttributeValue);
  REQUIRE(token.value == "1.0");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::DeclarationEnd);

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::TagBegin);
  REQUIRE(token.value == "Foo");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::AttributeName);
  REQUIRE(token.value == "test");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::Equals);

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::AttributeValue);
  REQUIRE(token.value == "1");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::TagEnd);

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::TagBegin);
  REQUIRE(token.value == "Bar");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::AttributeName);
  REQUIRE(token.value == "test");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::Equals);

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::AttributeValue);
  REQUIRE(token.value == "2");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::TagEnd);

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::Content);
  REQUIRE(token.value == "String\n");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::TagBegin);
  REQUIRE(token.value == "Baz");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::TagSelfClose);

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::TagBegin);
  REQUIRE(token.value == "Baz");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::AttributeName);
  REQUIRE(token.value == "attr");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::Equals);

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::AttributeValue);
  REQUIRE(token.value == "Hello World!");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::TagSelfClose);

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::TagClose);
  REQUIRE(token.value == "Bar");

  token = tokenizer.next_token().ok_value();
  REQUIRE(token.type == XML::TokenType::TagClose);
  REQUIRE(token.value == "Foo");

  REQUIRE(tokenizer.done());
}

TEST_CASE("Test XML parsing", "[IO]")
{
  std::stringstream xml;
  xml << "<?xml version=\"1.0\"?>\n";
  xml << "<Foo test=\"1\">\n";
  xml << "<Bar test=\"2\">\n";
  xml << "String\n";
  xml << "<Baz/>\n";
  xml << "<Baz attr = \"Hello World!\" />\n";
  xml << "</Bar  >\n";
  xml << "</Foo>\n";

  Result<std::unique_ptr<XML::XMLNode>, std::string> parse_result = XML::parse_xml(xml);

  REQUIRE(parse_result.is_ok());

  std::unique_ptr<XML::XMLNode> root = std::move(parse_result).take_ok();

  REQUIRE(root->attributes.size() == 1);
  REQUIRE(root->attributes.count("version") == 1);
  REQUIRE(root->attributes["version"] == "1.0");

  REQUIRE(root->children.size() == 1);
  XML::XMLNode* foo = root->children.front().get();

  REQUIRE(foo->name == "Foo");
  REQUIRE(foo->attributes.size() == 1);
  REQUIRE(foo->attributes.count("test") == 1);
  REQUIRE(foo->attributes["test"] == "1");

  REQUIRE(foo->children.size() == 1);
  XML::XMLNode* bar = foo->children.front().get();

  REQUIRE(bar->name == "Bar");
  REQUIRE(bar->attributes.size() == 1);
  REQUIRE(bar->attributes.count("test") == 1);
  REQUIRE(bar->attributes["test"] == "2");

  REQUIRE(bar->children.size() == 2);
  XML::XMLNode* baz_one = bar->children.front().get();
  XML::XMLNode* baz_two = bar->children.back().get();

  REQUIRE(baz_one->attributes.empty());
  REQUIRE(baz_two->attributes.size() == 1);
  REQUIRE(baz_two->attributes.count("attr") == 1);
  REQUIRE(baz_two->attributes["attr"] == "Hello World!");

  REQUIRE(baz_one->children.empty());
  REQUIRE(baz_two->children.empty());
}
