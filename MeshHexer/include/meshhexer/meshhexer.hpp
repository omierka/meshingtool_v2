#pragma once

#include <meshhexer/types.hpp>

#include <memory>

namespace MeshHexer
{
  /**
   * \brief MeshHexer surface mesh class
   *
   * Main API-entry-point for MeshHexer.
   * Contains a triangulated surface and can be used to either calculate properties of the surface mesh
   * or create volume meshes from the surface (non-aligned meshes only, for now).
   *
   * To create a SurfaceMesh instance, load a triangle mesh from disk using the MeshHexer::load_from_file() method.
   */
  class SurfaceMesh
  {
  public:
    // Forward declaration for PIMPL
    class SurfaceMeshImpl;

  private:

    /// Pointer to implementation
    std::unique_ptr<SurfaceMeshImpl> impl;

  public:
    /// Constructor
    explicit SurfaceMesh(std::unique_ptr<SurfaceMeshImpl> ptr);

    /// Destructor
    ~SurfaceMesh();

    /// Move-constructor
    SurfaceMesh(SurfaceMesh&&) noexcept;

    /// Copy-constructor
    SurfaceMesh(const SurfaceMesh&) = delete;

    /// Move-assign operator
    SurfaceMesh& operator=(SurfaceMesh&&) noexcept;

    /// Copy-assign operator
    SurfaceMesh& operator=(const SurfaceMesh&) = delete;

    /// standard ctor
    SurfaceMesh();


    /**
     * \brief Axis-aligned bounding box of the mesh
     *
     * \returns The axis-aligned bounding box of the mesh
     */
    BoundingBox bounding_box() const;

    /// Returns the number of vertices in the surface mesh
    std::uint32_t num_vertices() const;
    /// Returns the number of edges in the surface mesh
    std::uint32_t num_edges() const;
    /// Returns the number of faces in the surface mesh
    std::uint32_t num_faces() const;

    /**
     * \brief Check if the surface mesh describes a closed surface
     *
     * \returns True, if there are no border edges in the surface mesh
     */
    bool is_closed() const;

    /**
     * \brief Check if the winding order of the faces of the surface mesh is consistent
     *
     * \return True, if either all faces of the surface mesh are wound clockwise or all faces of the surface mesh are wound counter-clockwise
     */
    bool is_wound_consistently() const;

    /**
     * \brief Check if the surface normals of all faces of the surface mesh point towards the unbounded side
     *
     * \pre The surface mesh is closed, call SurfaceMesh::is_closed() to check
     * \pre The surface mesh is wound consistently, call SurfaceMesh::is_wound_consistently() to check
     *
     * \returns True, if the surface normals of all faces of the surface mesh point towards the unbounded side, i.e. outside.
     */
    bool is_outward_oriented() const;

    /**
     * \brief Compute the smallest aspect ratio of any face of the surface mesh
     *
     * \returns The minimal aspect ratio of any face of the surface mesh
     */
    double minimal_aspect_ratio() const;

    /**
     * \brief Compute the largest aspect ratio of any face of the surface mesh
     *
     * \returns The maximum aspect ratio of any face of the surface mesh
     */
    double maximal_aspect_ratio() const;

    /**
     * \brief Compute Gap candidates of the surface mesh.
     *
     * Gaps try to capture the diameter of the volume described by the surface mesh.
     * The diameter is measured from the barycentric centers of all faces of the surface mesh.
     * For each face of the surface mesh the largest sphere that just touches the center of
     * the face and one other point of the surface mesh is computed.
     * We call this sphere the maximal inscribed sphere (MIS) and use the diameter of that sphere
     * as the diameter at that face.
     *
     * \verbatim
     * │                 │
     * │     ∙∙∙∙∙∙∙     │
     * │  ∙∙∙       ∙∙∙  │
     * │ ∙             ∙ │
     * │∙    Diameter   ∙│
     * │∙---------------∙│
     * │∙               ∙│
     * │ ∙             ∙ │
     * │  ∙∙∙       ∙∙∙  │
     * │     ∙∙∙∙∙∙∙     │
     * │                 │
     * \endverbatim
     *
     * Note that not all these MIS span gaps that correspond to the intuitive diameter of the
     * surface mesh at that point. Badly reconstructed meshes for example might contain many
     * small "valleys" on surfaces that are supposed to be "smooth"
     *
     * \verbatim
     *              ∙∙∙∙∙∙∙
     *            ∙∙       ∙∙
     *\         ∙           ∙         /
     *  \\      ∙             ∙      //
     *    \\    ∙             ∙    //
     *      \\  ∙             ∙  //
     *        \\ ∙           ∙ //
     *          \\∙∙       ∙∙//
     *            \\∙∙∙∙∙∙∙//
     *              \\   //
     *                \ /
     * \endverbatim
     *
     * Each gap is given a confidence score in the range [0, 1] that shows how likely that gap
     * corresponds to a real diameter. That score is determined from:
     * - the aspect ratios of involved faces
     * - self-intersections of the surface mesh
     * - relative sizes of involved faces
     * - similarity of normals at the touching points of the MIS
     * - topological distance (distance across the surface between the touching points of the MIS)
     * - dihedral angles near the touching points of the MIS
     *
     * Generally gaps with a score above 0.95 have a good change of corresponding to real diameters.
     *
     * \returns A vector containing one gap candiate for each face of the surface mesh
     */
    std::vector<Gap> gaps();

    /**
     * \brief Compute the smallest Gap of the surface mesh
     *
     * See SurfaceMesh::gaps() for details on how gaps area measured and scored.
     *
     * \returns The smallest gap, by diameter, of all gaps with a score above 0.95,
     * or the smallest gap among the top 10 percent of all gaps, if no such gaps exist.
     */
    Gap min_gap();

    /**
     * \brief Create a non-fitted volume mesh for this surface mesh.
     *
     * \param[in] levels Number of regular refinements in intended multi-grid hierarchy
     * \param[in] settings Settings for the FBM mesh generator
     *
     * The volume mesh is created such that after the adaptive refinement prescribed
     * by the volume mesh and after \c settings.levels further global refinement steps,
     * the final volume mesh roughly hits the local min gap of the surface mesh.
     *
     * Assumes that adaptive refinements are 3-refinements and global refinements are 2-refinements.
     */
    VolumeMesh fbm_mesh(const FBMMeshSettings& settings);

    /**
     * \brief Create warnings for this surface mesh
     *
     * Warns about:
     * - self-intersections,
     * - degenerate triangles,
     * - highly anisotropic triangles,
     *
     * \returns Lists of warnings, separated by types
     */
    MeshWarnings warnings() const;

    /**
     * \brief Write the surface mesh to disk
     *
     * \param[in] filename Filename to write to. Must end in .ply or .vtu
     *
     * Writes the surface mesh to disk as a .ply or .vtu file. The written file
     * contains mesh properties that have been calculated as intermediate results,
     * such as maximal inscribed spheres, or topological distances.
     * These properties will be reused if the written mesh is read again.
     *
     * \returns A result indicating success or containing an error message
     */
    Result<void, std::string> write_to_file(const std::string& filename);
  };

  /**
   * \brief Load a surface mesh from a file
   *
   * \param[in] filename     Path to mesh file
   * \param[in] triangulate  If set, input mesh will automatically be triangulated
   *
   * Supported file formats are
   * - .off
   * - .obj
   * - .stl
   * - .ply
   * - .ts
   * - .vtp
   * - .vtu
   *
   * \returns a result containing either the surface mesh or an error message
   */
  Result<SurfaceMesh, std::string>
  load_from_file(const std::string& filename, bool triangulate = false);

} // namespace MeshHexer
