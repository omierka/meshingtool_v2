#pragma once

#include <cgal_types.hpp>
#include <meshhexer/types.hpp>

#include <initializer_list>
#include <queue>

namespace MeshHexer
{
  namespace XML
  {
    /**
     * \brief Enumeration of possible XML tokens
     */
    enum class TokenType
    {
      DeclarationBegin, // <?
      DeclarationEnd,   // ?>
      TagBegin,         // <
      TagEnd,           // >
      TagClose,         // </Tag>
      TagSelfClose,     // />
      AttributeName,    //
      Equals,           // =
      AttributeValue,
      Content,          // Anything between TagEnd / TagSelfClose and TagBegin
    };

    /// Output operator for TokenType
    inline std::ostream& operator<<(std::ostream& s, TokenType type)
    {
      switch(type)
      {
      case TokenType::DeclarationBegin: s << "TokenType::DeclarationBegin"; break;
      case TokenType::DeclarationEnd:   s << "TokenType::DeclarationEnd"; break;
      case TokenType::TagBegin:         s << "TokenType::TagBegin"; break;
      case TokenType::TagEnd:           s << "TokenType::TagEnd"; break;
      case TokenType::TagClose:         s << "TokenType::TagClose"; break;
      case TokenType::TagSelfClose:     s << "TokenType::TagSelfClose"; break;
      case TokenType::AttributeName:    s << "TokenType::AttributeName"; break;
      case TokenType::Equals:           s << "TokenType::AttEquals"; break;
      case TokenType::AttributeValue:   s << "TokenType::AttributeValue"; break;
      case TokenType::Content:          s << "TokenType::Content"; break;
      default:                          break;
      }
      return s;
    }

    /**
     * \brief XML-Token
     *
     * Contains the type of the token and optionally some value.
     * Value contains tag names, attribute names, attribute values, and
     * content of content tokens.
     */
    struct Token
    {
      TokenType type;
      std::string value;

      Token(TokenType t, std::string&& s) : type(t), value(std::move(s))
      {
      }
    };

    /**
     * \brief XML Tokenizer
     *
     * Converts an input stream containing an XML document into a series of XML tokens
     */
    class XMLTokenizer
    {
      /// Stream to read XML from
      std::istream& _stream;

      /// Current line read from stream
      std::string _line;

      /// Token buffer
      std::queue<Token> _tokens;

    public:
      /// Constructor
      explicit XMLTokenizer(std::istream& s);

      /**
       * \brief Get next token
       *
       * \returns a result containing either the next token in the XML document,
       * or a string describing an error.
       */
      Result<Token, std::string> next_token();

      /**
       * \brief Have all tokens been read
       *
       * \returns True, if the XMLTokenizer can produce no further tokens
       */
      bool done();

    private:
      /**
       * \brief Parse function
       *
       * Parses another segment of the XML document and pushes the generated
       * tokens into the token queue.
       */
      Result<void, std::string> parse();

      /**
       * \brief Parse a XML tag
       *
       * Handles parsing of XML tags and attributes, closing tags, xml-declarations, and comments.
       */
      Result<void, std::string> parse_tag();

      /**
       * \brief Parse XML attributes
       *
       * Parses attributes of an XML tag
       */
      Result<void, std::string> parse_attribute();

      /// Find the first occurence of a given char in the remaining XML document
      std::string::size_type find_first_of(char c);

      /// Find the first occurence of a given char in the first n chars of the remaining XML document
      std::optional<std::string::size_type> find_first_of_in_n(char c, std::size_t limit);

      /// Consume the XML document until the first occurence of any char in s
      std::string consume_until_any_of(const char* s);

      /// Consume n chars of the XML document
      std::string consume_n(std::string::size_type n);

      /// Consume one char of the XML document
      void consume_char();

      /// Look at the next char of the XML document
      char peek();
    };

    /**
     * \brief Class representing a node in a XML tree
     *
     * For simplicity, all attributes and any content are kept as strings:
     */
    struct XMLNode
    {
      /// Tag of this node
      std::string name;

      /// Attributes of this node
      std::map<std::string, std::string> attributes;

      /// Content of this node, if any
      std::vector<std::string> content;

      /// Children of this node, if any
      std::vector<std::unique_ptr<XMLNode>> children;

      /**
       * \brief Get children with given name
       *
       * \param[in] child_name Name of children
       *
       * \returns a vector of pointers to children with the given name
       */
      std::vector<XMLNode*> children_by_name(const std::string& child_name);

      /**
       * \brief Get child with given name
       *
       * \param[in] child_name Name of child
       *
       * \returns a pointer to the first child with the given name
       */
      XMLNode* child_by_name(const std::string& child_name);

      /**
       * \brief Get node at given path
       *
       * \param[in] path List of tag names
       *
       * The input list is treated as a path in the XML tree.
       * The first string of the path is looked up in the root element.
       * The second string of the path in that node, and so on
       */
      XMLNode* get_child(std::initializer_list<std::string> path);
    };

    /**
     * \brief Parse an XML document into a tree of XML nodes
     *
     * \param[in] stream Stream containing a XML document
     * \returns a result containing either the root of the XML-tree or an error string
     */
    Result<std::unique_ptr<XMLNode>, std::string> parse_xml(std::istream& stream);
  } // namespace XML

  void write_geo_compound_2d(const std::string& filename, const Polygon2D& poly);
  void write_geo(const std::string& filename, const Polygon2D& poly);
  void write_polygon_brep(const std::string& filename, const Polygon2D& poly);
  void write_polygon(const std::string& filename, const Polygon2D& poly);
  void write_polygon(const std::string& filename, const PolygonWithHoles2D& poly);
  void write_polylines(const std::string& filename, const Polylines3D& polylines);
  void write_polylines(const std::string& filename, const Polylines2D& polylines);

  void write_feat_xml(std::ostream& stream, const VolumeMesh& vmesh);

  void write_vtu(std::ostream& stream, const Mesh& mesh);
  Result<Mesh, std::string> read_vtu(std::istream& stream);
} // namespace MeshHexer
