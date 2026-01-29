#include <cctype>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <stack>
#include <string>
#include <string_view>
#include <type_traits>

#include <cgal_types.hpp>
#include <io.hpp>
#include <macros.hpp>
#include <meshhexer/types.hpp>

namespace MeshHexer
{
  namespace
  {
    template<typename T>
    struct VTUTypeTrait
    {
      static constexpr std::string_view as_string = "";
    };

    template<>
    struct VTUTypeTrait<std::int8_t>
    {
      static constexpr std::string_view as_string = "Int8";
    };

    template<>
    struct VTUTypeTrait<std::int16_t>
    {
      static constexpr std::string_view as_string = "Int16";
    };

    template<>
    struct VTUTypeTrait<std::int32_t>
    {
      static constexpr std::string_view as_string = "Int32";
    };

    template<>
    struct VTUTypeTrait<std::int64_t>
    {
      static constexpr std::string_view as_string = "Int64";
    };

    template<>
    struct VTUTypeTrait<std::uint8_t>
    {
      static constexpr std::string_view as_string = "UInt8";
    };

    template<>
    struct VTUTypeTrait<std::uint16_t>
    {
      static constexpr std::string_view as_string = "UInt16";
    };

    template<>
    struct VTUTypeTrait<std::uint32_t>
    {
      static constexpr std::string_view as_string = "UInt32";
    };

    template<>
    struct VTUTypeTrait<std::uint64_t>
    {
      static constexpr std::string_view as_string = "UInt64";
    };

    template<>
    struct VTUTypeTrait<float>
    {
      static constexpr std::string_view as_string = "Float32";
    };

    template<>
    struct VTUTypeTrait<double>
    {
      static constexpr std::string_view as_string = "Float64";
    };

    template<typename I, typename T>
    void vtu_write_property(std::ostream& stream, const Mesh& mesh, const std::string& prop_name)
    {
      using PropertyMap = typename Mesh::Property_map<I, T>;

      std::optional<PropertyMap> mmap = mesh.property_map<I, T>(prop_name);

      if(mmap.has_value())
      {
        PropertyMap map = mmap.value();

        std::string_view type = VTUTypeTrait<T>::as_string;

        stream << "<DataArray type=\"" << type << "\" Name=\"" << prop_name.substr(2) << "\" format=\"ascii\">\n";

        for(const T& val : map)
        {
          stream << val << "\n";
        }

        stream << "</DataArray>\n";
      }
    }

    template<typename I>
    void vtu_write_property(std::ostream& stream, const Mesh& mesh, const std::string& prop_name)
    {
      vtu_write_property<I, std::int8_t>(stream, mesh, prop_name);
      vtu_write_property<I, std::int16_t>(stream, mesh, prop_name);
      vtu_write_property<I, std::int32_t>(stream, mesh, prop_name);
      vtu_write_property<I, std::int64_t>(stream, mesh, prop_name);
      vtu_write_property<I, std::uint8_t>(stream, mesh, prop_name);
      vtu_write_property<I, std::uint16_t>(stream, mesh, prop_name);
      vtu_write_property<I, std::uint32_t>(stream, mesh, prop_name);
      vtu_write_property<I, std::uint64_t>(stream, mesh, prop_name);
      vtu_write_property<I, float>(stream, mesh, prop_name);
      vtu_write_property<I, double>(stream, mesh, prop_name);
    }

    void trim(std::string& s)
    {
      const auto not_whitespace = [](char c) { return !std::isspace(c); };

      // Erase from beginning to first non-whitespace character
      s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_whitespace));

      // Erase from last non-whitespace character to end
      s.erase(std::find_if(s.rbegin(), s.rend(), not_whitespace).base(), s.end());
    }

    template<typename T>
    Result<std::vector<T>, std::string> parse_content(const XML::XMLNode* n)
    {
      std::vector<T> values;
      for(const std::string& content : n->content)
      {
        const char* s = content.c_str();
        const char* end = s + content.size();
        while(*s != 0)
        {
          while(std::isspace(*s))
          {
            s++;
          }

          if(*s == 0)
          {
            break;
          }

          T value;
          std::from_chars_result result = std::from_chars(s, end, value);

          if(result.ec == std::errc::invalid_argument)
          {
            return Result<std::vector<T>, std::string>::err("Conversion failed: invalid argument");
          }
          if(result.ec == std::errc::result_out_of_range)
          {
            return Result<std::vector<T>, std::string>::err("Conversion failed: result out of range");
          }

          values.push_back(value);
          s = result.ptr;
        }
      }

      return Result<std::vector<T>, std::string>::ok(std::move(values));
    }

    template<typename T>
    Result<void, std::string> parse_mesh_property(Mesh& m, const XML::XMLNode* n)
    {
      std::string name = std::string("f:") + n->attributes.at("Name");
      auto map = m.add_property_map<FaceIndex, T>(name).first;

      Result<std::vector<T>, std::string> value_result = parse_content<T>(n);

      if(value_result.is_err())
      {
        return Result<void, std::string>(std::move(value_result).take_err());
      }

      std::vector<T> vals = std::move(value_result).take_ok();

      for(std::size_t i(0); i < vals.size(); i++)
      {
        map[FaceIndex(i)] = vals[i];
      }

      return Result<void, std::string>::ok();
    }
  } // namespace

  namespace XML
  {
    XMLTokenizer::XMLTokenizer(std::istream& s) : _stream(s)
    {
      std::getline(_stream, _line);
      trim(_line);
    }

    /**
     * \brief Get next token of XML document
     */
    Result<Token, std::string> XMLTokenizer::next_token()
    {
      using ResultType = Result<Token, std::string>;

      while(_tokens.empty() && (!_line.empty() || _stream.good()))
      {
        Result<void, std::string> parse_result = parse();
        if(parse_result.is_err())
        {
          return ResultType::err(parse_result.err_value());
        }
      }

      Token t = _tokens.front();
      _tokens.pop();

      return ResultType::ok(t);
    }

    Result<void, std::string> XMLTokenizer::parse()
    {
      if(peek() == '<')
      {
        parse_tag();
      }
      else
      {
        _tokens.emplace(TokenType::Content, consume_until_any_of("<"));
      }

      return Result<void, std::string>::ok();
    }

    Result<void, std::string> XMLTokenizer::parse_tag()
    {
      // Consume leading '<'
      consume_char();

      char next = peek();
      if(next == '/')
      {
        // Found start of a closing tag

        // Consume the '/'
        consume_char();

        std::string::size_type pos = find_first_of('>');
        if(pos == std::string::npos)
        {
          return Result<void, std::string>::err("Expected closing '>' for closing tag");
        }

        _tokens.emplace(TokenType::TagClose, consume_n(pos));

        XASSERT(peek() == '>');
        consume_char();
      }
      else if(next == '?')
      {
        // Consume ?xml
        consume_char();
        consume_char();
        consume_char();
        consume_char();

        _tokens.emplace(TokenType::DeclarationBegin, "");

        // Found start of a declaration tag
        Result<void, std::string> attribute_result = parse_attribute();

        if(attribute_result.is_err())
        {
          return attribute_result;
        }

        std::string::size_type pos = find_first_of('>');
        std::string suffix = consume_n(pos + 1);

        if(suffix != "?>")
        {
          return Result<void, std::string>::err("Expected '?>' after attributes of XML-declaration");
        }

        _tokens.emplace(TokenType::DeclarationEnd, "");
      }
      else if(next == '!')
      {
        // Consume !--
        consume_char();
        consume_char();
        consume_char();

        while(true)
        {
          std::string::size_type pos = find_first_of('-');
          consume_n(pos);

          if(peek() == '-')
          {
            consume_char();
            if(peek() == '-')
            {
              consume_char();
              if(peek() == '>')
              {
                consume_char();
                break;
              }
            }
          }
        }
      }
      else
      {
        // Found start of a begin tag

        // Find tag name
        std::string tag = consume_until_any_of(" \t\n\r\f\v/>");
        trim(tag);

        _tokens.emplace(TokenType::TagBegin, std::move(tag));

        // Parse attributes, if any
        Result<void, std::string> attribute_result = parse_attribute();
        if(attribute_result.is_err())
        {
          return attribute_result;
        }

        if(peek() == '/')
        {
          _tokens.emplace(TokenType::TagSelfClose, "");

          // Consume "/>"
          consume_n(2);
        }
        else if(peek() == '>')
        {
          _tokens.emplace(TokenType::TagEnd, "");

          // Consume '>'
          consume_char();
        }
      }

      return Result<void, std::string>::ok();
    }

    Result<void, std::string> XMLTokenizer::parse_attribute()
    {
      std::string::size_type end_of_tag = find_first_of('>');
      std::string::size_type next_equals = find_first_of_in_n('=', end_of_tag).value_or(std::string::npos);

      while(next_equals != std::string::npos && next_equals < end_of_tag)
      {
        std::string attribute_name = consume_n(next_equals);
        trim(attribute_name);

        ASSERT(peek() == '=');
        consume_char();
        ASSERT(peek() == '"');
        consume_char();

        std::string::size_type end_quote = find_first_of('"');

        std::string attribute_value = consume_n(end_quote);
        trim(attribute_value);

        _tokens.emplace(TokenType::AttributeName, std::move(attribute_name));
        _tokens.emplace(TokenType::Equals, "");
        _tokens.emplace(TokenType::AttributeValue, std::move(attribute_value));

        ASSERT(peek() == '"');
        consume_char();

        end_of_tag = find_first_of('>');
        next_equals = find_first_of_in_n('=', end_of_tag).value_or(std::string::npos);
      }

      return Result<void, std::string>::ok();
    }

    char XMLTokenizer::peek()
    {
      while(_line.empty() && _stream.good())
      {
        std::getline(_stream, _line);
        trim(_line);
      }
      trim(_line);

      XASSERT(_line.empty() || !std::isspace(_line.front()));

      return _line.front();
    }

    std::string::size_type XMLTokenizer::find_first_of(char c)
    {
      std::string::size_type result = 0;
      const char* s = _line.c_str();

      std::string tmp;
      while(!_line.empty() || _stream.good())
      {
        while(*s != 0)
        {
          if(*s == c)
          {
            return result;
          }
          result++;
          s++;
        }

        if(_stream.good())
        {
          std::getline(_stream, tmp);
          _line.append(tmp + "\n");
        }

        // Continue searching where we left off
        s = _line.c_str() + result;
      }
      return result;
    }

    std::optional<std::string::size_type> XMLTokenizer::find_first_of_in_n(char c, std::size_t limit)
    {
      std::string::size_type result = 0;

      const char* s = _line.c_str();

      while(true)
      {
        while(*s != 0 && (limit == 0 || result < limit))
        {
          if(*s == c)
          {
            return result;
          }
          result++;
          s++;
        }

        if(result == limit)
        {
          return std::nullopt;
        }

        if(_stream.good())
        {
          std::string tmp;
          std::getline(_stream, tmp);
          _line.append(tmp + "\n");
        }
        else
        {
          break;
        }

        // Continue searching where we left off
        s = _line.c_str() + result;
      }
      return result;
    }

    std::string XMLTokenizer::consume_n(std::string::size_type n)
    {
      std::vector<std::string> lines;
      std::size_t size = 0;

      while(size < n && (!_line.empty() || _stream.good()))
      {
        if(_line.size() < n - size)
        {
          // Consume whole line
          size += _line.size();
          lines.push_back(_line);

          if(_stream.good())
          {
            std::getline(_stream, _line);
            _line.append("\n");
          }
        }
        else
        {
          // Take part of line
          std::string::size_type diff = n - size;
          lines.push_back(_line.substr(0, diff));
          _line = _line.substr(diff);
          size += diff;
        }
      }

      trim(_line);

      XASSERT(_line.empty() || !std::isspace(_line.front()));

      return std::accumulate(lines.begin(), lines.end(), std::string{});
    }

    std::string XMLTokenizer::consume_until_any_of(const char* s)
    {
      std::vector<std::string> lines;

      while(!_line.empty() || _stream.good())
      {
        std::string::size_type pos = _line.find_first_of(s);

        if(pos == std::string::npos)
        {
          // Consume whole line
          lines.push_back(_line);
          _line.clear();

          if(_stream.good())
          {
            std::getline(_stream, _line);
            _line.append("\n");
          }
        }
        else
        {
          // Take part of line and stop consuming
          lines.push_back(_line.substr(0, pos));
          _line = _line.substr(pos);
          break;
        }
      }

      trim(_line);

      XASSERT(_line.empty() || !std::isspace(_line.front()));

      return std::accumulate(lines.begin(), lines.end(), std::string{});
    }

    void XMLTokenizer::consume_char()
    {
      while(_line.empty() && _stream.good())
      {
        std::string tmp;
        std::getline(_stream, tmp);
        _line.append(tmp + "\n");
      }
      _line = _line.substr(1);
      trim(_line);

      while(_line.empty() && _stream.good())
      {
        std::getline(_stream, _line);
      }

      XASSERT(_line.empty() || !std::isspace(_line.front()));
    }

    bool XMLTokenizer::done()
    {
      return _line.empty() && !_stream.good() && _tokens.empty();
    }

    Result<std::unique_ptr<XMLNode>, std::string> parse_xml(std::istream& stream)
    {
      using ResultType = Result<std::unique_ptr<XMLNode>, std::string>;

      XMLTokenizer tokenizer(stream);

      auto root = std::make_unique<XMLNode>();

      std::stack<XMLNode*> node_stack;
      node_stack.push(root.get());

      std::string attribute_name;

      while(!tokenizer.done())
      {
        Result<Token, std::string> token_result = tokenizer.next_token();

        if(token_result.is_err())
        {
          return ResultType::err(std::move(token_result).take_err());
        }

        Token t = token_result.ok_value();

        switch(t.type)
        {
        case TokenType::TagBegin:
          {
            auto node = std::make_unique<XMLNode>();
            node->name = t.value;

            XMLNode* parent = node_stack.top();
            node_stack.push(node.get());
            parent->children.push_back(std::move(node));
            break;
          }
        case TokenType::AttributeName:
          {
            attribute_name = t.value;
            break;
          }
        case TokenType::AttributeValue:
          {
            node_stack.top()->attributes[attribute_name] = t.value;
            break;
          }
        case TokenType::Content:
          {
            node_stack.top()->content.push_back(t.value);
            break;
          }
        case TokenType::DeclarationBegin:
          {
            if(node_stack.size() != 1)
            {
              return ResultType::err("XML-declaration must be at root-scope!");
            }
            break;
          }
        case TokenType::TagClose:
        case TokenType::TagSelfClose:
          {
            node_stack.pop();
            break;
          }
        case TokenType::TagEnd:
        case TokenType::Equals:
        case TokenType::DeclarationEnd:
          {
            // Nothing to do
            break;
          }
        }
      }

      return ResultType::ok(std::move(root));
    }

    std::vector<XMLNode*> XMLNode::children_by_name(const std::string& child_name)
    {
      std::vector<XMLNode*> result;
      for(const auto& child : children)
      {
        if(child->name == child_name)
        {
          result.push_back(child.get());
        }
      }

      return result;
    }

    XMLNode* XMLNode::child_by_name(const std::string& child_name)
    {
      for(const auto& child : children)
      {
        if(child->name == child_name)
        {
          return child.get();
        }
      }
      return nullptr;
    }

    XMLNode* XMLNode::get_child(std::initializer_list<std::string> path)
    {
      XMLNode* result = this;
      for(const std::string& s : path)
      {
        result = result->child_by_name(s);

        if(result == nullptr)
        {
          return nullptr;
        }
      }

      return result;
    }
  } // namespace XML

  void write_geo_compound_2d(const std::string& filename, const Polygon2D& poly)
  {
    std::ofstream output(filename);

    if(!output)
    {
      std::cerr << "Could not open file " << filename << " for writing.\n";
      return;
    }

    // Set options
    output << "Mesh.Algorithm = 8;\n";              // Frontal-Delaunay for Quads
    output << "Mesh.RecombinationAlgorithm = 3;\n"; // Blossom Full-Quad
    output << "Mesh.Format = 16;\n";                // vtk
    output << "Mesh.RecombineAll = 1;\n";           // Always produce quad meshes

    int next_tag = 1;

    for(const Point2D& p : poly)
    {
      output << "Point(" << next_tag++ << ") = {" << p.x() << ", " << p.y() << ", 0, 5.0};\n";
    }

    next_tag = 1;
    for(std::size_t i(0); i < poly.size() - 1; i++)
    {
      output << "Curve(" << next_tag++ << ") = {" << std::to_string(i + 1) << ", " << std::to_string(i + 2) << "};\n";
    }
    output << "Curve(" << next_tag++ << ") = {" << poly.size() << ", 1};\n";

    int line_loop_tag = next_tag;
    output << "Curve Loop(" << next_tag++ << ") = {";
    for(std::size_t i(0); i < poly.size(); i++)
    {
      output << std::to_string(i + 1);
      if(i + 1 != poly.size())
      {
        output << ", ";
      }
    }
    output << "};\n";

    int surface_tag = next_tag;
    output << "Plane Surface(" << next_tag++ << ") = {" << line_loop_tag << "};\n";

    output << "Compound Curve{";
    for(std::size_t i(0); i < poly.size(); i++)
    {
      output << std::to_string(i + 1);
      if(i + 1 != poly.size())
      {
        output << ", ";
      }
    }
    output << "};\n";
    output << "Compound Surface{" << surface_tag << "};\n";
    output << "Mesh 2\n";
  }

  void write_geo(const std::string& filename, const Polygon2D& poly)
  {
    std::ofstream output(filename);

    if(!output)
    {
      std::cerr << "Could not open file " << filename << " for writing.\n";
      return;
    }

    // Set options
    output << "Mesh.Algorithm = 8;\n";              // Frontal-Delaunay for Quads
    output << "Mesh.RecombinationAlgorithm = 3;\n"; // Blossom Full-Quad
    output << "Mesh.Format = 16;\n";                // vtk
    output << "Mesh.RecombineAll = 1;\n";           // Always produce quad meshes

    std::size_t next_tag = 1;

    for(const Point2D& p : poly)
    {
      output << "Point(" << next_tag++ << ") = {" << p.x() << ", " << p.y() << ", 0, 5.0};\n";
    }

    next_tag = 1;
    for(std::size_t i(0); i < poly.size() - 1; i++)
    {
      output << "Line(" << next_tag++ << ") = {" << std::to_string(i + 1) << ", " << std::to_string(i + 2) << "};\n";
    }
    output << "Line(" << next_tag++ << ") = {" << poly.size() << ", 1};\n";

    std::size_t line_loop_tag = next_tag;
    output << "Line Loop(" << next_tag++ << ") = {";
    for(std::size_t i(0); i < poly.size(); i++)
    {
      output << std::to_string(i + 1);
      if(i + 1 != poly.size())
      {
        output << ", ";
      }
    }
    output << "};\n";

    // std::size_t surface_tag = next_tag;
    output << "Plane Surface(" << next_tag++ << ") = {" << line_loop_tag << "};\n";

    // output << "Extrude {0, 0, 10} {Surface{" << surface_tag << "}; Layers
    // {1}; }\n";
    output << "Mesh 2\n";
  }

  void write_polygon_brep(const std::string& filename, const Polygon2D& poly)
  {
    std::ofstream output(filename);

    if(!output)
    {
      std::cerr << "Could not open file " << filename << " for writing.\n";
      return;
    }

    output << "DBRep_DrawableShape\n\n";
    output << "CASCADE Topology V1, (c)  Matra-Datavision\n";
    output << "Locations 0\n";
    output << "Curve2ds 0\n";
    output << "Curves 0\n";
    output << "Polygon3D 1\n";                       // One 3D polygon
    output << std::to_string(poly.size()) << " 0\n"; // Number of points and parameter presence. 0 = no parameters.
    output << "0.1\n";                               // deflection

    // Points of the polygon, all on a single line, in triplets
    for(const Point2D& p : poly)
    {
      output << p.x() << " " << p.y() << " 0 ";
    }
    output << "\n";

    output << "PolygonOnTriangulations 0\n";
    output << "Surfaces 0\n";
    output << "Triangulations 0\n";
    output << "TShapes 0\n";
  }

  void write_polygon(const std::string& filename, const Polygon2D& poly)
  {
    std::ofstream output(filename);

    if(output)
    {
      output << "<VTKFile type=\"PolyData\">\n";
      output << "<PolyData>\n";
      output << "<Piece NumberOfPoints=\"" << poly.size()
             << "\" NumberOfVerts=\"0\" NumberOfLines=\"1\" "
                "NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
      output << "<PointData></PointData>\n";
      output << "<CellData></CellData>\n";
      output << "<Points>\n";
      output << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
                "Format=\"ascii\">\n";

      for(const Point2D& p : poly.vertices())
      {
        output << p.x() << " " << p.y() << " 0\n";
      }

      output << "</DataArray>\n";
      output << "</Points>\n";
      output << "<Verts></Verts>\n";
      output << "<Lines>\n";
      output << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";

      for(std::size_t i(0); i < poly.size(); i++)
      {
        output << unsigned(i) << " ";
      }
      output << "0\n";
      output << "</DataArray>\n";

      output << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
      output << poly.size() + 1 << "\n";
      output << "</DataArray>\n";

      output << "</Lines>\n";
      output << "<Strips></Strips>\n";
      output << "<Polys></Polys>\n";
      output << "</Piece>\n";
      output << "</PolyData>\n";
      output << "</VTKFile>\n";
    }
    else
    {
      std::cerr << "Failed to write polygon to " << filename << "!\n";
    }
  }

  void write_polygon(const std::string& filename, const PolygonWithHoles2D& poly)
  {
    std::ofstream output(filename);

    if(output)
    {
      std::size_t num_polys = 1 + poly.holes().size();
      std::size_t total_points = 0;
      std::vector<std::size_t> starting_points;

      starting_points.push_back(total_points);
      total_points += poly.outer_boundary().size();
      for(const Polygon2D& hole : poly.holes())
      {
        starting_points.push_back(total_points);
        total_points += hole.size();
      }

      output << "<VTKFile type=\"PolyData\">\n";
      output << "<PolyData>\n";
      output << "<Piece NumberOfPoints=\"" << total_points << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << num_polys
             << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
      output << "<PointData></PointData>\n";
      output << "<CellData></CellData>\n";
      output << "<Points>\n";
      output << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
                "Format=\"ascii\">\n";

      for(const Point2D& p : poly.outer_boundary().vertices())
      {
        output << p.x() << " " << p.y() << " 0\n";
      }

      for(const Polygon2D& hole : poly.holes())
      {
        for(const Point2D& p : hole.vertices())
        {
          output << p.x() << " " << p.y() << " 0\n";
        }
      }

      output << "</DataArray>\n";
      output << "</Points>\n";
      output << "<Verts></Verts>\n";
      output << "<Lines>\n";
      output << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";

      if(poly.outer_boundary().size() > 0)
      {
        int idx = 0;
        for(std::size_t i(0); i < poly.outer_boundary().size() - 1; i++)
        {
          output << unsigned(i) << " ";
        }
        output << starting_points[idx] << "\n";
        idx++;

        for(const Polygon2D& hole : poly.holes())
        {
          std::size_t starting_index = starting_points[idx];
          for(std::size_t i(0); i < hole.size() - 1; i++)
          {
            output << starting_index + i << "\n";
          }
          output << starting_index << "\n";
          idx++;
        }
      }

      output << "</DataArray>\n";

      output << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
      for(std::size_t i(1); i < starting_points.size(); i++)
      {
        output << starting_points[i] << "\n";
      }
      output << total_points << "\n";
      output << "</DataArray>\n";

      output << "</Lines>\n";
      output << "<Strips></Strips>\n";
      output << "<Polys></Polys>\n";
      output << "</Piece>\n";
      output << "</PolyData>\n";
      output << "</VTKFile>\n";
    }
    else
    {
      std::cerr << "Failed to write polygon to " << filename << "!\n";
    }
  }

  void write_polylines(const std::string& filename, const Polylines3D& polylines)
  {
    std::ofstream output(filename);

    if(output)
    {
      std::size_t num_lines = polylines.size();
      std::size_t total_points = 0;
      std::vector<std::size_t> starting_points;

      for(const Polyline3D& line : polylines)
      {
        starting_points.push_back(total_points);
        total_points += line.size();
      }

      output << "<VTKFile type=\"PolyData\">\n";
      output << "<PolyData>\n";
      output << "<Piece NumberOfPoints=\"" << total_points << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << num_lines
             << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
      output << "<PointData></PointData>\n";
      output << "<CellData></CellData>\n";
      output << "<Points>\n";
      output << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
                "Format=\"ascii\">\n";

      for(const Polyline3D& polyline : polylines)
      {
        for(const Point3D& p : polyline)
        {
          output << p.x() << " " << p.y() << " " << p.z() << "\n";
        }
      }

      output << "</DataArray>\n";
      output << "</Points>\n";
      output << "<Verts></Verts>\n";
      output << "<Lines>\n";
      output << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";

      auto it = polylines.begin();
      for(std::size_t i(0); i < polylines.size(); i++)
      {
        std::size_t starting_idx = starting_points[i];

        for(std::size_t j(0); j < it->size() - 1; j++)
        {
          output << starting_idx + j << "\n";
        }
        output << starting_idx << "\n";
        it++;
      }
      output << "</DataArray>\n";

      output << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
      for(std::size_t i(1); i < starting_points.size(); i++)
      {
        output << starting_points[i] << "\n";
      }
      output << total_points << "\n";
      output << "</DataArray>\n";

      output << "</Lines>\n";
      output << "<Strips></Strips>\n";
      output << "<Polys></Polys>\n";
      output << "</Piece>\n";
      output << "</PolyData>\n";
      output << "</VTKFile>\n";
    }
    else
    {
      std::cerr << "Failed to write polyline to " << filename << "!\n";
    }
  }

  void write_polylines(const std::string& filename, const Polylines2D& polylines)
  {
    std::ofstream output(filename);

    if(output)
    {
      std::size_t num_lines = polylines.size();
      std::size_t total_points = 0;
      std::vector<std::size_t> starting_points;

      for(const Polyline2D& line : polylines)
      {
        starting_points.push_back(total_points);
        total_points += line.size();
      }

      output << "<VTKFile type=\"PolyData\">\n";
      output << "<PolyData>\n";
      output << "<Piece NumberOfPoints=\"" << total_points << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << num_lines
             << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
      output << "<PointData></PointData>\n";
      output << "<CellData></CellData>\n";
      output << "<Points>\n";
      output << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
                "Format=\"ascii\">\n";

      for(const Polyline2D& polyline : polylines)
      {
        for(const Point2D& p : polyline)
        {
          output << p.x() << " " << p.y() << " 0\n";
        }
      }

      output << "</DataArray>\n";
      output << "</Points>\n";
      output << "<Verts></Verts>\n";
      output << "<Lines>\n";
      output << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";

      auto it = polylines.begin();
      for(std::size_t i(0); i < polylines.size(); i++)
      {
        std::size_t starting_idx = starting_points[i];

        for(std::size_t j(0); j < it->size() - 1; j++)
        {
          output << starting_idx + j << "\n";
        }
        output << starting_idx << "\n";
        it++;
      }
      output << "</DataArray>\n";

      output << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
      for(std::size_t i(1); i < starting_points.size(); i++)
      {
        output << starting_points[i] << "\n";
      }
      output << total_points << "\n";
      output << "</DataArray>\n";

      output << "</Lines>\n";
      output << "<Strips></Strips>\n";
      output << "<Polys></Polys>\n";
      output << "</Piece>\n";
      output << "</PolyData>\n";
      output << "</VTKFile>\n";
    }
    else
    {
      std::cerr << "Failed to write polyline to " << filename << "!\n";
    }
  }

  void write_vtu(std::ostream& stream, const Mesh& mesh)
  {
    const std::size_t num_verts = mesh.number_of_vertices();
    const std::size_t num_faces = mesh.number_of_faces();

    stream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    stream << "<UnstructuredGrid>\n";
    stream << "<Piece NumberOfPoints=\"" << num_verts << "\" NumberOfCells=\"" << num_faces << "\">\n";
    stream << "<Points>\n";
    stream << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for(const VertexIndex v : mesh.vertices())
    {
      stream << mesh.point(v) << "\n";
    }
    stream << "</DataArray>\n";
    stream << "</Points>\n";
    stream << "<Cells>\n";
    stream << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";

    for(const FaceIndex f : mesh.faces())
    {
      for(const VertexIndex v : mesh.vertices_around_face(mesh.halfedge(f)))
      {
        stream << static_cast<std::int64_t>(v) << "\n";
      }
    }
    stream << "</DataArray>\n";
    stream << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";

    for(std::size_t i(0); i < num_faces; i++)
    {
      stream << (i + 1) * 3 << "\n";
    }

    stream << "</DataArray>\n";
    stream << "<DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";

    constexpr int triangle_type = 5;
    for(std::size_t i(0); i < num_faces; i++)
    {
      stream << triangle_type << "\n";
    }
    stream << "</DataArray>\n";
    stream << "</Cells>\n";
    stream << "<CellData>\n";

    for(const std::string& property_name : mesh.properties<VertexIndex>())
    {
      vtu_write_property<VertexIndex>(stream, mesh, property_name);
    }
    for(const std::string& property_name : mesh.properties<HalfedgeIndex>())
    {
      vtu_write_property<HalfedgeIndex>(stream, mesh, property_name);
    }
    for(const std::string& property_name : mesh.properties<EdgeIndex>())
    {
      vtu_write_property<EdgeIndex>(stream, mesh, property_name);
    }
    for(const std::string& property_name : mesh.properties<FaceIndex>())
    {
      vtu_write_property<FaceIndex>(stream, mesh, property_name);
    }

    stream << "</CellData>\n";
    stream << "</Piece>\n";
    stream << "</UnstructuredGrid>\n";
    stream << "</VTKFile>\n";
  }

  Result<Mesh, std::string> read_vtu(std::istream& stream)
  {
    using ResultType = Result<Mesh, std::string>;

    Result<std::unique_ptr<XML::XMLNode>, std::string> parse_result = XML::parse_xml(stream);
    if(parse_result.is_err())
    {
      Result<Mesh, std::string>::err(std::move(parse_result).take_err());
    }

    std::unique_ptr<XML::XMLNode> root = std::move(parse_result).take_ok();
    Mesh result;

    XML::XMLNode* piece_node = root->get_child({"VTKFile", "UnstructuredGrid", "Piece"});
    if(piece_node == nullptr)
    {
      return ResultType::err("Could not find node with path \"VTKFile/UnstructuredGrid/Piece\"");
    }

    XML::XMLNode* points_node = piece_node->get_child({"Points", "DataArray"});
    if(points_node == nullptr)
    {
      return ResultType::err("Could not find node with path \"VTKFile/UnstructuredGrid/Piece/Points/DataArray\"");
    }

    if(points_node->attributes.count("NumberOfComponents") == 0)
    {
      return ResultType::err(
        "Expected attribute 'NumberOfComponents' at node 'VTKFile/UnstructuredGrid/Piece/Points/DataArray'");
    }

    if(points_node->attributes["NumberOfComponents"] != "3")
    {
      return ResultType::err("Only 3D meshes are supported");
    }

    auto coords_result = parse_content<double>(points_node);
    if(coords_result.is_err())
    {
      return ResultType::err(std::move(coords_result).take_err());
    }
    std::vector<double> coords = std::move(coords_result).take_ok();

    ASSERT(coords.size() % 3 == 0);

    std::vector<VertexIndex> vertices;
    for(std::size_t i(0); i < coords.size(); i += 3)
    {
      vertices.push_back(result.add_vertex(Point3D{coords[i], coords[i + 1], coords[i + 2]}));
    }

    // Parse topology

    XML::XMLNode* cell_node = piece_node->child_by_name("Cells");

    if(cell_node == nullptr)
    {
      return ResultType::err("Could not find node with path 'VTKFile/UnstructuredGrid/Piece/Cells'");
    }

    std::vector<XML::XMLNode*> array_nodes = cell_node->children_by_name("DataArray");
    std::vector<const XML::XMLNode*> property_nodes;

    const XML::XMLNode* connectivity_node = nullptr;
    const XML::XMLNode* offsets_node = nullptr;
    const XML::XMLNode* types_node = nullptr;

    for(const XML::XMLNode* n : array_nodes)
    {
      if(n->attributes.at("Name") == "connectivity")
      {
        connectivity_node = n;
      }
      else if(n->attributes.at("Name") == "offsets")
      {
        offsets_node = n;
      }
      else if(n->attributes.at("Name") == "types")
      {
        types_node = n;
      }
    }

    if(connectivity_node == nullptr)
    {
      return ResultType::err("Could not find mesh connectivity");
    }
    if(offsets_node == nullptr)
    {
      return ResultType::err("Could not find mesh connectivity offsets");
    }
    if(types_node == nullptr)
    {
      return ResultType::err("Could not find mesh element types");
    }

    auto connectivity_result = parse_content<std::int64_t>(connectivity_node);
    if(connectivity_result.is_err())
    {
      return ResultType::err(std::move(connectivity_result).take_err());
    }
    std::vector<std::int64_t> connectivity = std::move(connectivity_result).take_ok();

    auto offset_result = parse_content<std::int64_t>(offsets_node);
    if(offset_result.is_err())
    {
      return ResultType::err(std::move(offset_result).take_err());
    }
    std::vector<std::int64_t> offsets = std::move(offset_result).take_ok();

    auto type_result = parse_content<std::int64_t>(types_node);
    if(type_result.is_err())
    {
      return ResultType::err(std::move(type_result).take_err());
    }
    std::vector<std::int64_t> types = std::move(type_result).take_ok();

    for(std::int64_t type : types)
    {
      XASSERT(type == 5);
    }

    for(std::int64_t i : offsets)
    {
      auto offset = static_cast<std::size_t>(i);

      std::size_t a = connectivity[offset - 3];
      std::size_t b = connectivity[offset - 2];
      std::size_t c = connectivity[offset - 1];

      FaceIndex f = result.add_face(vertices[a], vertices[b], vertices[c]);

      if(f == Mesh::null_face())
      {
        return ResultType::err("Failed to create face");
      }
    }

    // Parse properties

    XML::XMLNode* cell_data_node = piece_node->child_by_name("CellData");

    if(cell_data_node != nullptr)
    {
      for(const XML::XMLNode* n : cell_data_node->children_by_name("DataArray"))
      {
        std::string type = n->attributes.at("type");

        if(type == "UInt8")
        {
          parse_mesh_property<std::uint8_t>(result, n);
        }
        else if(type == "UInt16")
        {
          parse_mesh_property<std::uint16_t>(result, n);
        }
        else if(type == "UInt32")
        {
          parse_mesh_property<std::uint32_t>(result, n);
        }
        else if(type == "UInt64")
        {
          parse_mesh_property<std::uint64_t>(result, n);
        }
        else if(type == "Int8")
        {
          parse_mesh_property<std::int8_t>(result, n);
        }
        else if(type == "Int16")
        {
          parse_mesh_property<std::int16_t>(result, n);
        }
        else if(type == "Int32")
        {
          parse_mesh_property<std::int32_t>(result, n);
        }
        else if(type == "Int64")
        {
          parse_mesh_property<std::int64_t>(result, n);
        }
        else if(type == "Float32")
        {
          parse_mesh_property<float>(result, n);
        }
        else if(type == "Float64")
        {
          parse_mesh_property<double>(result, n);
        }
        else
        {
          return ResultType::err(std::string("Unknown type for mesh property: " + type));
        }
      }
    }

    return Result<Mesh, std::string>::ok(std::move(result));
  }
} // namespace MeshHexer
