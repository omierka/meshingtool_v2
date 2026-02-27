module def_mod
  use iso_fortran_env, only: real64
  use iso_c_binding, only: c_ptr, c_null_ptr, c_char, c_null_char, c_size_t, c_double, c_int, c_associated
  use mpi
  use setupe3dfile_reader, only: MeshConfig
  use preprocessor_config_mod, only: get_tolerance_floor, get_meshdeform_characteristic_scale
  implicit none
  private

  integer, parameter, public :: rk = real64
  integer, parameter, public :: constraint_none = 0
  integer, parameter, public :: constraint_x_min = 1
  integer, parameter, public :: constraint_x_max = 2
  integer, parameter, public :: constraint_y_min = 4
  integer, parameter, public :: constraint_y_max = 8
  integer, parameter, public :: constraint_z_min = 16
  integer, parameter, public :: constraint_z_max = 32
  integer, parameter, public :: constraint_cyl_inner = 64
  integer, parameter, public :: constraint_cyl_outer = 128

  character(len=*), parameter :: constraint_names(8) = [ &
       character(len=16) :: 'x-', 'x+', 'y-', 'y+', 'z-', 'z+', 'inner', 'outer' &
  ]
  integer, parameter :: constraint_masks(8) = [ &
       constraint_x_min, constraint_x_max, constraint_y_min, constraint_y_max, &
       constraint_z_min, constraint_z_max, constraint_cyl_inner, constraint_cyl_outer &
  ]
  integer, parameter :: hexahedron_edges(2, 12) = reshape([ &
       1, 2, 2, 3, 3, 4, 4, 1, 5, 6, 6, 7, 7, 8, 8, 5, 1, 5, 2, 6, 3, 7, 4, 8 &
     ], [2, 12])
  integer, parameter :: hexahedron_faces(4, 6) = reshape([ &
       1, 2, 3, 4, &
       5, 6, 7, 8, &
       1, 2, 6, 5, &
       2, 3, 7, 6, &
       3, 4, 8, 7, &
       4, 1, 5, 8], [4, 6])

  type, public :: hex_mesh_type
     integer :: nel = 0
     integer :: nvt = 0
     integer :: nbct = 0
     integer :: nve = 0
     integer :: nee = 0
     integer :: nae = 0
     integer :: nedge = 0
     integer :: edge_capacity = 0
     real(rk), allocatable :: dcorvg(:,:)
     integer, allocatable :: kvert(:,:)
     integer, allocatable :: knpr(:)
     integer, allocatable :: kedge(:,:)
     integer, allocatable :: edge_nodes(:,:)
   contains
     procedure :: clear => clear_hex_mesh
  end type hex_mesh_type

  type, public :: deformation_context
     logical :: initialized = .false.
  end type deformation_context

  type, public :: vertex_constraint_lists
     integer, allocatable :: x_min(:)
     integer, allocatable :: x_max(:)
     integer, allocatable :: y_min(:)
     integer, allocatable :: y_max(:)
     integer, allocatable :: z_min(:)
     integer, allocatable :: z_max(:)
     integer, allocatable :: cyl_inner(:)
     integer, allocatable :: cyl_outer(:)
     integer, allocatable :: boundary(:)
  end type vertex_constraint_lists

  type, public :: surface_mesh
     type(c_ptr) :: handle = c_null_ptr
     integer :: vertex_count = 0
     integer :: triangle_count = 0
     logical :: aabb_ready = .false.
     character(len=1024) :: source_file = ''
   contains
     procedure :: clear => clear_surface_mesh
  end type surface_mesh

  public :: initialize_deformation_context
  public :: finalize_deformation_context
  public :: summarize_mesh
  public :: assign_vertex_constraints
  public :: report_constraint_summary
  public :: build_constraint_lists
  public :: release_constraint_lists
  public :: report_constraint_combinations
  public :: identify_boundary_nodes
  public :: build_edge_to_element
  public :: broadcast_mesh
  public :: broadcast_mesh_coordinates
  public :: report_boundary_assignment_status
  public :: load_surface_mesh
  public :: report_surface_summary
  public :: apply_edge_deformation
  public :: compute_signed_distances
  public :: compute_signed_distances_parallel
  public :: report_distance_summary
  public :: write_deformed_vtu
  public :: write_deformed_tri

  interface
     function cgal_load_off(path) bind(C, name="cgal_load_off") result(handle)
       import :: c_ptr, c_char
       character(kind=c_char), dimension(*) :: path
       type(c_ptr) :: handle
     end function cgal_load_off

     subroutine cgal_free_mesh(handle) bind(C, name="cgal_free_mesh")
       import :: c_ptr
       type(c_ptr), value :: handle
     end subroutine cgal_free_mesh

     function cgal_get_vertex_count(handle) bind(C, name="cgal_get_vertex_count") result(count)
       import :: c_ptr, c_size_t
       type(c_ptr), value :: handle
       integer(c_size_t) :: count
     end function cgal_get_vertex_count

     function cgal_get_triangle_count(handle) bind(C, name="cgal_get_triangle_count") result(count)
       import :: c_ptr, c_size_t
       type(c_ptr), value :: handle
       integer(c_size_t) :: count
     end function cgal_get_triangle_count

     function cgal_build_aabb_tree(handle) bind(C, name="cgal_build_aabb_tree") result(status)
       import :: c_ptr, c_int
       type(c_ptr), value :: handle
       integer(c_int) :: status
     end function cgal_build_aabb_tree

     function cgal_signed_distance_to_mesh(handle, point, distance_out) bind(C, name="cgal_signed_distance_to_mesh") result(status)
       import :: c_ptr, c_double, c_int
       type(c_ptr), value :: handle
       real(c_double), intent(in) :: point(*)
       real(c_double), intent(out) :: distance_out
       integer(c_int) :: status
     end function cgal_signed_distance_to_mesh
  end interface

contains

  subroutine initialize_deformation_context(context)
    type(deformation_context), intent(inout) :: context

    context%initialized = .true.
  end subroutine initialize_deformation_context

  subroutine finalize_deformation_context(context, mesh)
    type(deformation_context), intent(inout) :: context
    type(hex_mesh_type), intent(inout), optional :: mesh

    if (present(mesh)) call mesh%clear()
    context%initialized = .false.
  end subroutine finalize_deformation_context

  subroutine clear_hex_mesh(mesh)
    class(hex_mesh_type), intent(inout) :: mesh

    mesh%nel = 0
    mesh%nvt = 0
    mesh%nbct = 0
    mesh%nve = 0
    mesh%nee = 0
    mesh%nae = 0
    mesh%nedge = 0
    mesh%edge_capacity = 0
    if (allocated(mesh%dcorvg)) deallocate(mesh%dcorvg)
    if (allocated(mesh%kvert))  deallocate(mesh%kvert)
    if (allocated(mesh%knpr))   deallocate(mesh%knpr)
    if (allocated(mesh%kedge))  deallocate(mesh%kedge)
    if (allocated(mesh%edge_nodes)) deallocate(mesh%edge_nodes)
  end subroutine clear_hex_mesh

  subroutine clear_surface_mesh(surface)
    class(surface_mesh), intent(inout) :: surface

    if (c_associated(surface%handle)) call cgal_free_mesh(surface%handle)
    surface%handle = c_null_ptr
    surface%vertex_count = 0
    surface%triangle_count = 0
    surface%aabb_ready = .false.
    surface%source_file = ''
  end subroutine clear_surface_mesh

  subroutine summarize_mesh(mesh, label)
    type(hex_mesh_type), intent(in) :: mesh
    character(len=*), intent(in), optional :: label
    real(rk) :: bounds_min(3), bounds_max(3)
    integer :: ivt

    if (mesh%nel <= 0 .or. mesh%nvt <= 0) then
       write(*, '(A)') 'Mesh contains no elements; nothing to report.'
       return
    end if

    bounds_min = huge(1.0_rk)
    bounds_max = -huge(1.0_rk)
    if (allocated(mesh%dcorvg)) then
       do ivt = 1, size(mesh%dcorvg, 2)
          bounds_min = min(bounds_min, mesh%dcorvg(:, ivt))
          bounds_max = max(bounds_max, mesh%dcorvg(:, ivt))
       end do
    end if

    if (present(label)) then
       write(*, '(A)') 'Mesh summary for: ' // trim(label)
    else
       write(*, '(A)') 'Mesh summary:'
    end if
    write(*, '(A,I0)') '  Number of elements (NEL): ', mesh%nel
    write(*, '(A,I0)') '  Number of vertices (NVT): ', mesh%nvt
    write(*, '(A,I0)') '  Node boundary codes (NBCT): ', mesh%nbct
    write(*, '(A,I0)') '  Vertices per element (NVE): ', mesh%nve
    write(*, '(A,I0)') '  Element edges (NEE): ', mesh%nee
    write(*, '(A,I0)') '  Element faces (NAE): ', mesh%nae
    if (allocated(mesh%dcorvg)) then
       write(*, '(A,3(1X,ES13.6))') '  Bounding box min: ', bounds_min
       write(*, '(A,3(1X,ES13.6))') '  Bounding box max: ', bounds_max
    end if
  end subroutine summarize_mesh

  subroutine assign_vertex_constraints(mesh, config, constraint_flags)
    type(hex_mesh_type), intent(in) :: mesh
    type(MeshConfig), intent(in) :: config
    integer, allocatable, intent(out) :: constraint_flags(:)
    character(len=:), allocatable :: mesh_type

    if (mesh%nvt <= 0) then
       if (allocated(constraint_flags)) deallocate(constraint_flags)
       allocate(constraint_flags(0))
       return
    end if

    if (.not.allocated(mesh%dcorvg)) then
       if (allocated(constraint_flags)) deallocate(constraint_flags)
       allocate(constraint_flags(mesh%nvt))
       constraint_flags = constraint_none
       return
    end if

    if (allocated(constraint_flags)) deallocate(constraint_flags)
    allocate(constraint_flags(mesh%nvt))
    constraint_flags = constraint_none

    mesh_type = lowercase(trim(config%mesh_type))
    select case (trim(mesh_type))
    case ('box')
       call assign_box_constraints(mesh, config, constraint_flags, compute_min_edge_length(mesh))
    case ('hollowcylinder')
       call assign_cylinder_constraints(mesh, config, constraint_flags, compute_min_edge_length(mesh))
    case default
       ! Unknown mesh type: keep all vertices unconstrained.
    end select
  end subroutine assign_vertex_constraints

  subroutine assign_box_constraints(mesh, config, constraint_flags, min_edge_length)
    type(hex_mesh_type), intent(in) :: mesh
    type(MeshConfig), intent(in) :: config
    integer, intent(inout) :: constraint_flags(:)
    real(rk), intent(in) :: min_edge_length
    real(rk) :: start(3), finish(3), lengths(3)
    real(rk) :: tol, coord(3), max_dim, base_tol
    integer :: ivt

    call compute_mesh_bounds(mesh, start, finish)
    if (config%box%has_geometry_start) start = real(config%box%geometry_start, rk)
    if (config%box%has_geometry_length) finish = start + real(config%box%geometry_length, rk)

    lengths = max(finish - start, 0.0_rk)
    max_dim = max(1.0_rk, maxval(lengths))
    if (min_edge_length > 0.0_rk) then
       tol = max(real(get_tolerance_floor(), rk), 0.25_rk * min_edge_length)
    else
       tol = max(real(get_tolerance_floor(), rk), 1.0e-6_rk * max_dim)
    end if

    do ivt = 1, mesh%nvt
       coord = mesh%dcorvg(:, ivt)
       if (abs(coord(1) - start(1)) <= tol) constraint_flags(ivt) = ior(constraint_flags(ivt), constraint_x_min)
       if (abs(coord(1) - finish(1)) <= tol) constraint_flags(ivt) = ior(constraint_flags(ivt), constraint_x_max)
       if (abs(coord(2) - start(2)) <= tol) constraint_flags(ivt) = ior(constraint_flags(ivt), constraint_y_min)
       if (abs(coord(2) - finish(2)) <= tol) constraint_flags(ivt) = ior(constraint_flags(ivt), constraint_y_max)
       if (abs(coord(3) - start(3)) <= tol) constraint_flags(ivt) = ior(constraint_flags(ivt), constraint_z_min)
       if (abs(coord(3) - finish(3)) <= tol) constraint_flags(ivt) = ior(constraint_flags(ivt), constraint_z_max)
    end do
  end subroutine assign_box_constraints

  subroutine assign_cylinder_constraints(mesh, config, constraint_flags, min_edge_length)
    type(hex_mesh_type), intent(in) :: mesh
    type(MeshConfig), intent(in) :: config
    integer, intent(inout) :: constraint_flags(:)
    real(rk), intent(in) :: min_edge_length
    real(rk) :: bounds_min(3), bounds_max(3)
    real(rk) :: center(2), coord(3)
    real(rk) :: axial_min, axial_max
    real(rk) :: inner_radius, outer_radius
    real(rk) :: tol_rad, tol_axial, base_tol
    logical :: have_inner, have_outer, have_axial_min, have_axial_max
    integer :: ivt

    call compute_mesh_bounds(mesh, bounds_min, bounds_max)
    center(1) = 0.5_rk * (bounds_min(1) + bounds_max(1))
    center(2) = 0.5_rk * (bounds_min(2) + bounds_max(2))

    have_outer = config%cylinder%has_barrel_diameter
    if (have_outer) then
       outer_radius = 0.5_rk * real(config%cylinder%barrel_diameter, rk)
    else
       outer_radius = 0.0_rk
    end if
    have_inner = config%cylinder%has_inner_diameter
    if (have_inner) then
       inner_radius = 0.5_rk * real(config%cylinder%inner_diameter, rk)
    else
       inner_radius = 0.0_rk
    end if
    have_axial_min = config%cylinder%has_axial_start
    if (have_axial_min) then
       axial_min = real(config%cylinder%axial_start, rk)
    else
       axial_min = bounds_min(3)
    end if
    have_axial_max = have_axial_min .and. config%cylinder%has_barrel_length
    if (have_axial_max) then
       axial_max = real(config%cylinder%axial_start + config%cylinder%barrel_length, rk)
    else
       axial_max = bounds_max(3)
    end if

    if (min_edge_length > 0.0_rk) then
       tol_rad = max(real(get_tolerance_floor(), rk), 0.25_rk * min_edge_length)
       tol_axial = tol_rad
    else
       base_tol = max(real(get_tolerance_floor(), rk), 1.0e-6_rk)
       tol_rad = max(base_tol, 1.0e-6_rk * max(max(outer_radius, inner_radius), 1.0_rk))
       tol_axial = max(base_tol, 1.0e-6_rk * max(axial_max - axial_min, 1.0_rk))
    end if

    do ivt = 1, mesh%nvt
       coord = mesh%dcorvg(:, ivt)
       if (have_outer) then
          if (abs(radial_distance(coord, center) - outer_radius) <= tol_rad) then
             constraint_flags(ivt) = ior(constraint_flags(ivt), constraint_cyl_outer)
          end if
       end if
       if (have_inner) then
          if (abs(radial_distance(coord, center) - inner_radius) <= tol_rad) then
             constraint_flags(ivt) = ior(constraint_flags(ivt), constraint_cyl_inner)
          end if
       end if
       if (have_axial_min) then
          if (abs(coord(3) - axial_min) <= tol_axial) &
               constraint_flags(ivt) = ior(constraint_flags(ivt), constraint_z_min)
       end if
       if (have_axial_max) then
          if (abs(coord(3) - axial_max) <= tol_axial) &
               constraint_flags(ivt) = ior(constraint_flags(ivt), constraint_z_max)
       end if
    end do
  end subroutine assign_cylinder_constraints

  subroutine load_surface_mesh(surface, filename)
    type(surface_mesh), intent(inout) :: surface
    character(len=*), intent(in) :: filename
    character(len=:), allocatable :: trimmed
    character(kind=c_char), allocatable :: c_path(:)
    logical :: exists
    integer(c_int) :: status

    trimmed = trim(filename)
    inquire(file=trimmed, exist=exists)
    if (.not.exists) then
       write(*, '(A)') 'Surface mesh file not found: ' // trimmed
       stop 1
    end if

    call surface%clear()
    c_path = to_c_string(trimmed)
    surface%handle = cgal_load_off(c_path)
    if (.not.c_associated(surface%handle)) then
       write(*, '(A)') 'Failed to load surface mesh from: ' // trimmed
       stop 1
    end if

    surface%vertex_count = int(cgal_get_vertex_count(surface%handle))
    surface%triangle_count = int(cgal_get_triangle_count(surface%handle))
    status = cgal_build_aabb_tree(surface%handle)
    if (status /= 0) then
       write(*, '(A,I0)') 'Failed to prepare surface AABB tree. STATUS=', status
       stop 1
    end if
    surface%aabb_ready = .true.
    surface%source_file = trimmed
    if (allocated(c_path)) deallocate(c_path)
  end subroutine load_surface_mesh

  subroutine report_surface_summary(surface)
    type(surface_mesh), intent(in) :: surface

    if (.not.c_associated(surface%handle)) then
      write(*, '(A)') 'Surface mesh not loaded.'
      return
    end if

    write(*, '(A)') 'Surface mesh summary:'
    if (len_trim(surface%source_file) > 0) then
       write(*, '(A)') '  Source: ' // trim(surface%source_file)
    end if
    write(*, '(A,I0)') '  Vertices: ', surface%vertex_count
    write(*, '(A,I0)') '  Triangles: ', surface%triangle_count
    if (surface%aabb_ready) then
       write(*, '(A)') '  AABB tree: ready'
    else
       write(*, '(A)') '  AABB tree: not ready'
    end if
  end subroutine report_surface_summary

  subroutine report_constraint_summary(constraint_flags, boundary_mask)
    integer, intent(in) :: constraint_flags(:)
    logical, intent(in), optional :: boundary_mask(:)
    integer :: total, boundary_total, unassigned

    total = size(constraint_flags)
    if (total <= 0) then
       write(*, '(A)') 'No vertices available for constraint summary.'
       return
    end if

    call print_constraint_count('x-', constraint_x_min, constraint_flags)
    call print_constraint_count('x+', constraint_x_max, constraint_flags)
    call print_constraint_count('y-', constraint_y_min, constraint_flags)
    call print_constraint_count('y+', constraint_y_max, constraint_flags)
    call print_constraint_count('z-', constraint_z_min, constraint_flags)
    call print_constraint_count('z+', constraint_z_max, constraint_flags)
    call print_constraint_count('inner-cylinder', constraint_cyl_inner, constraint_flags)
    call print_constraint_count('outer-cylinder', constraint_cyl_outer, constraint_flags)
    call report_constraint_combinations(constraint_flags)

    if (present(boundary_mask)) then
       boundary_total = count(boundary_mask)
       write(*, '(A,I0)') 'Topological boundary vertices: ', boundary_total
       unassigned = count(boundary_mask .and. (constraint_flags == constraint_none))
       if (unassigned > 0) then
          write(*, '(A,I0)') 'WARNING: boundary vertices without parametrization: ', unassigned
       end if
    end if
  end subroutine report_constraint_summary

  subroutine print_constraint_count(label, bitmask, constraint_flags)
    character(len=*), intent(in) :: label
    integer, intent(in) :: bitmask
    integer, intent(in) :: constraint_flags(:)
    integer :: count_value

    count_value = count(iand(constraint_flags, bitmask) /= 0)
    if (count_value > 0) then
       write(*, '(A,1X,I0)') trim(label) // ' vertices:', count_value
    end if
  end subroutine print_constraint_count

  subroutine compute_mesh_bounds(mesh, bounds_min, bounds_max)
    type(hex_mesh_type), intent(in) :: mesh
    real(rk), intent(out) :: bounds_min(3)
    real(rk), intent(out) :: bounds_max(3)
    integer :: ivt

    bounds_min = huge(1.0_rk)
    bounds_max = -huge(1.0_rk)
    if (.not.allocated(mesh%dcorvg)) return
    do ivt = 1, size(mesh%dcorvg, 2)
       bounds_min = min(bounds_min, mesh%dcorvg(:, ivt))
       bounds_max = max(bounds_max, mesh%dcorvg(:, ivt))
    end do
  end subroutine compute_mesh_bounds

  real(rk) function radial_distance(coord, center) result(dist)
    real(rk), intent(in) :: coord(3)
    real(rk), intent(in) :: center(2)
    real(rk) :: dx, dy
    dx = coord(1) - center(1)
    dy = coord(2) - center(2)
    dist = sqrt(dx * dx + dy * dy)
  end function radial_distance

  pure function lowercase(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code
    do i = 1, len(text)
       code = iachar(text(i:i))
       if (code >= iachar('A') .and. code <= iachar('Z')) then
          out(i:i) = achar(code + 32)
       else
          out(i:i) = text(i:i)
       end if
    end do
  end function lowercase

  function to_c_string(str) result(c_chars)
    character(len=*), intent(in) :: str
    character(kind=c_char), allocatable :: c_chars(:)
    integer :: n, i

    n = len_trim(str)
    allocate(c_chars(n + 1))
    do i = 1, n
       c_chars(i) = str(i:i)
    end do
    c_chars(n + 1) = c_null_char
  end function to_c_string

  subroutine compute_signed_distances(surface, mesh, distances)
    type(surface_mesh), intent(in) :: surface
    type(hex_mesh_type), intent(in) :: mesh
    real(rk), allocatable, intent(out) :: distances(:)
    integer :: ivt
    real(c_double) :: point(3), signed_distance
    integer(c_int) :: status

    if (.not.c_associated(surface%handle)) then
       write(*, '(A)') 'Surface mesh handle is not available for distance queries.'
       allocate(distances(0))
       return
    end if

    if (mesh%nvt <= 0) then
       allocate(distances(0))
       return
    end if

    allocate(distances(mesh%nvt))
    do ivt = 1, mesh%nvt
       point(1) = real(mesh%dcorvg(1, ivt), c_double)
       point(2) = real(mesh%dcorvg(2, ivt), c_double)
       point(3) = real(mesh%dcorvg(3, ivt), c_double)
       status = cgal_signed_distance_to_mesh(surface%handle, point, signed_distance)
       if (status /= 0) then
          distances(ivt) = 0.0_rk
       else
          distances(ivt) = real(signed_distance, rk)
       end if
    end do
  end subroutine compute_signed_distances

  subroutine compute_signed_distances_parallel(mesh, surface, distances, mpi_rank, mpi_size)
    type(hex_mesh_type), intent(in) :: mesh
    type(surface_mesh), intent(in) :: surface
    real(rk), allocatable, intent(out) :: distances(:)
    integer, intent(in) :: mpi_rank, mpi_size

    integer :: nvt, start_idx, end_idx, local_count, idx, ivt
    real(rk), allocatable :: local_dist(:)
    integer, allocatable :: counts(:), displs(:)
    real(c_double) :: point(3), signed_distance
    integer :: mpi_err
    integer(c_int) :: status

    nvt = mesh%nvt
    if (nvt <= 0) then
       allocate(distances(0))
       return
    end if

    start_idx = (mpi_rank * nvt) / mpi_size + 1
    end_idx   = ((mpi_rank + 1) * nvt) / mpi_size
    local_count = max(0, end_idx - start_idx + 1)
    if (local_count > 0) then
       allocate(local_dist(local_count))
       do idx = 1, local_count
          ivt = start_idx + idx - 1
          point(1) = real(mesh%dcorvg(1, ivt), c_double)
          point(2) = real(mesh%dcorvg(2, ivt), c_double)
          point(3) = real(mesh%dcorvg(3, ivt), c_double)
          status = cgal_signed_distance_to_mesh(surface%handle, point, signed_distance)
          if (status /= 0) signed_distance = 0.0_c_double
          local_dist(idx) = real(signed_distance, rk)
       end do
    else
       allocate(local_dist(0))
    end if

    allocate(counts(mpi_size))
    allocate(displs(mpi_size))
    do idx = 0, mpi_size - 1
       counts(idx + 1) = ((idx + 1) * nvt) / mpi_size - (idx * nvt) / mpi_size
       displs(idx + 1) = ((idx * nvt) / mpi_size)
    end do

    allocate(distances(nvt))
    call MPI_Allgatherv(local_dist, local_count, MPI_DOUBLE_PRECISION, distances, counts, displs, &
         MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpi_err)

    if (allocated(local_dist)) deallocate(local_dist)
    if (allocated(counts)) deallocate(counts)
    if (allocated(displs)) deallocate(displs)
  end subroutine compute_signed_distances_parallel

  subroutine report_distance_summary(surface, mesh)
    type(surface_mesh), intent(in) :: surface
    type(hex_mesh_type), intent(in) :: mesh
    real(rk), allocatable :: distances(:)
    real(rk) :: min_val, max_val

    call compute_signed_distances(surface, mesh, distances)
    if (size(distances) <= 0) then
       write(*, '(A)') 'No signed distances computed.'
    else
       min_val = minval(distances)
       max_val = maxval(distances)
       write(*, '(A)') 'Signed distance summary:'
       write(*, '(A,1X,ES13.6)') '  minimum:', min_val
       write(*, '(A,1X,ES13.6)') '  maximum:', max_val
    end if
    if (allocated(distances)) deallocate(distances)
  end subroutine report_distance_summary

  subroutine apply_edge_deformation(mesh, config, constraint_flags, distances)
    type(hex_mesh_type), intent(inout) :: mesh
    type(MeshConfig), intent(in) :: config
    integer, intent(in) :: constraint_flags(:)
    real(rk), intent(in) :: distances(:)

    integer, parameter :: edge_vertices(2, 12) = reshape([ &
         1, 2, 2, 3, 3, 4, 4, 1, 1, 5, 2, 6, 3, 7, 4, 8, 5, 6, 6, 7, 7, 8, 8, 5], [2, 12])
    integer :: elem_idx, edge_idx, ivt1, ivt2
    integer :: nel, nvt
    real(rk) :: characteristic_size, scale_factor, weight_edge
    real(rk), allocatable :: volumes(:), vertex_measure(:), weights(:)
    real(rk), allocatable :: accum_x(:), accum_y(:), accum_z(:), accum_w(:), vertex_weights(:)
    real(rk) :: p1(3), p2(3), daux1, daux2, px, py, pz
    real(rk), parameter :: damping = 0.2_rk
    real(rk) :: box_start(3), box_finish(3)
    real(rk) :: cyl_center(2), inner_radius, outer_radius
    real(rk) :: axial_min, axial_max
    logical :: has_inner, has_outer, has_axial_min, has_axial_max

    nvt = mesh%nvt
    nel = mesh%nel
    if (nvt <= 0 .or. nel <= 0) return
    if (size(distances) /= nvt) then
       write(*, '(A)') 'Distance array size mismatch; deformation skipped.'
       return
    end if
    if (.not.allocated(mesh%kvert)) return
    if (.not.allocated(mesh%dcorvg)) return
    if (.not.allocated(mesh%kedge)) call build_edge_to_element(mesh)

    call compute_box_limits(mesh, config, box_start, box_finish)
    call compute_cylinder_parameters(mesh, config, cyl_center, inner_radius, outer_radius, axial_min, axial_max, &
         has_inner, has_outer, has_axial_min, has_axial_max)
    call compute_characteristic_size(mesh, config, characteristic_size)
    if (characteristic_size <= 0.0_rk) characteristic_size = 1.0_rk
    scale_factor = real(get_meshdeform_characteristic_scale(), rk) / characteristic_size

    allocate(volumes(nel))
    allocate(vertex_measure(nvt))
    allocate(weights(nvt))
    allocate(accum_x(nvt))
    allocate(accum_y(nvt))
    allocate(accum_z(nvt))
    allocate(accum_w(nvt))
    allocate(vertex_weights(nvt))

    call compute_vertex_weights(distances, scale_factor, vertex_weights)
    call compute_hexahedron_volumes(mesh, volumes)

    vertex_measure = 0.0_rk
    weights = 0.0_rk
    do elem_idx = 1, nel
       do edge_idx = 1, 8
          ivt1 = mesh%kvert(edge_idx, elem_idx)
          if (ivt1 < 1 .or. ivt1 > nvt) cycle
          vertex_measure(ivt1) = vertex_measure(ivt1) + abs(volumes(elem_idx))
          weights(ivt1) = weights(ivt1) + 1.0_rk
       end do
    end do
    do ivt1 = 1, nvt
       if (weights(ivt1) > 0.0_rk) vertex_measure(ivt1) = vertex_measure(ivt1) / weights(ivt1)
    end do

    do ivt1 = 1, nvt
       vertex_measure(ivt1) = vertex_measure(ivt1) * vertex_weights(ivt1)
    end do

    accum_x = 0.0_rk
    accum_y = 0.0_rk
    accum_z = 0.0_rk
    accum_w = 0.0_rk

    do elem_idx = 1, nel
       do edge_idx = 1, 12
          ivt1 = mesh%kvert(edge_vertices(1, edge_idx), elem_idx)
          ivt2 = mesh%kvert(edge_vertices(2, edge_idx), elem_idx)
          if (ivt1 < 1 .or. ivt1 > nvt) cycle
          if (ivt2 < 1 .or. ivt2 > nvt) cycle
          p1 = mesh%dcorvg(:, ivt1)
          p2 = mesh%dcorvg(:, ivt2)
          daux1 = abs(vertex_measure(ivt1))
          daux2 = abs(vertex_measure(ivt2))
          weight_edge = 1.0_rk

          accum_x(ivt1) = accum_x(ivt1) + weight_edge * p2(1) * daux2
          accum_y(ivt1) = accum_y(ivt1) + weight_edge * p2(2) * daux2
          accum_z(ivt1) = accum_z(ivt1) + weight_edge * p2(3) * daux2
          accum_w(ivt1) = accum_w(ivt1) + weight_edge * daux2

          accum_x(ivt2) = accum_x(ivt2) + weight_edge * p1(1) * daux1
          accum_y(ivt2) = accum_y(ivt2) + weight_edge * p1(2) * daux1
          accum_z(ivt2) = accum_z(ivt2) + weight_edge * p1(3) * daux1
          accum_w(ivt2) = accum_w(ivt2) + weight_edge * daux1
       end do
    end do

    do ivt1 = 1, nvt
       if (accum_w(ivt1) > 0.0_rk) then
          px = accum_x(ivt1) / accum_w(ivt1)
          py = accum_y(ivt1) / accum_w(ivt1)
          pz = accum_z(ivt1) / accum_w(ivt1)
          call enforce_parametrization_constraints(constraint_flags, ivt1, mesh%dcorvg(:, ivt1), px, py, pz, &
               box_start, box_finish, cyl_center, inner_radius, outer_radius, axial_min, axial_max, has_inner, &
               has_outer, has_axial_min, has_axial_max)
          mesh%dcorvg(1, ivt1) = max(0.0_rk, 1.0_rk - damping) * mesh%dcorvg(1, ivt1) + damping * px
          mesh%dcorvg(2, ivt1) = max(0.0_rk, 1.0_rk - damping) * mesh%dcorvg(2, ivt1) + damping * py
          mesh%dcorvg(3, ivt1) = max(0.0_rk, 1.0_rk - damping) * mesh%dcorvg(3, ivt1) + damping * pz
       end if
    end do

    deallocate(volumes, vertex_measure, weights, accum_x, accum_y, accum_z, accum_w, vertex_weights)
  end subroutine apply_edge_deformation

  subroutine write_deformed_vtu(mesh, surface, filename)
    type(hex_mesh_type), intent(in) :: mesh
    type(surface_mesh), intent(in) :: surface
    character(len=*), intent(in) :: filename
    real(rk), allocatable :: distances(:)

    call compute_signed_distances(surface, mesh, distances)
    call write_mesh_vtu(mesh, distances, filename)
    if (allocated(distances)) deallocate(distances)
  end subroutine write_deformed_vtu

  subroutine write_deformed_tri(mesh, filename)
    type(hex_mesh_type), intent(in) :: mesh
    character(len=*), intent(in) :: filename
    integer :: unit, nel, nvt, nbct, nve, nee, nae, i

    if (.not.allocated(mesh%dcorvg) .or. .not.allocated(mesh%kvert)) then
       write(*, '(A)') 'Mesh data unavailable; TRI file not written.'
       return
    end if

    nel = mesh%nel
    nvt = mesh%nvt
    nbct = mesh%nbct
    nve = mesh%nve
    nee = mesh%nee
    nae = mesh%nae

    open(newunit=unit, file=trim(filename), status='replace', action='write')
    write(unit, '(A)') 'MeshDeform output'
    write(unit, '(A)') 'Generated by meshdeform'
    write(unit, '(I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,3X,A)') nel, nvt, nbct, nve, nee, nae, 'NEL,NVT,NBCT,NVE,NEE,NAE'
    write(unit, '(A)') 'DCORVG'
    do i = 1, nvt
       write(unit, '(3(1X,ES24.16))') mesh%dcorvg(1, i), mesh%dcorvg(2, i), mesh%dcorvg(3, i)
    end do
    write(unit, '(A)') 'KVERT'
    do i = 1, nel
       write(unit, '(*(1X,I10))') mesh%kvert(:, i)
    end do
    write(unit, '(A)') 'KNPR'
    if (allocated(mesh%knpr)) then
       do i = 1, nvt
          write(unit, '(I8)') mesh%knpr(i)
       end do
    else
       do i = 1, nvt
          write(unit, '(I8)') 0
       end do
    end if
    close(unit)
    write(*, '(A,A)') 'TRI mesh written to ', trim(filename)
  end subroutine write_deformed_tri

  subroutine write_mesh_vtu(mesh, distances, filename)
    type(hex_mesh_type), intent(in) :: mesh
    real(rk), intent(in) :: distances(:)
    character(len=*), intent(in) :: filename

    integer :: unit, nvt, nel, nve, offset, i, j
    logical :: has_knpr

    if (.not.allocated(mesh%dcorvg) .or. .not.allocated(mesh%kvert)) then
       write(*, '(A)') 'Mesh data unavailable; VTU not written.'
       return
    end if

    nvt = size(mesh%dcorvg, 2)
    nel = size(mesh%kvert, 2)
    nve = size(mesh%kvert, 1)
    has_knpr = allocated(mesh%knpr)

    if (size(distances) /= nvt) then
       write(*, '(A)') 'Distance data size mismatch; VTU not written.'
       return
    end if

    open(newunit=unit, file=trim(filename), status='replace', action='write')
    write(unit, '(A)') '<?xml version="1.0"?>'
    write(unit, '(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(unit, '(A)') '  <UnstructuredGrid>'
    write(unit, '(A,I0,A,I0,A)') '    <Piece NumberOfPoints="', nvt, '" NumberOfCells="', nel, '">'

    write(unit, '(A)') '      <PointData>'
    if (has_knpr) then
       write(unit, '(A)') '        <DataArray type="Int32" Name="KNPR" format="ascii">'
       write(unit, '(*(1X,I8))') (int(mesh%knpr(i), kind=4), i = 1, nvt)
       write(unit, '(A)') '        </DataArray>'
    end if
    write(unit, '(A)') '        <DataArray type="Float64" Name="distance" format="ascii">'
    write(unit, '(*(1X,ES24.16))') (distances(i), i = 1, nvt)
    write(unit, '(A)') '        </DataArray>'
    write(unit, '(A)') '      </PointData>'

    write(unit, '(A)') '      <CellData/>'

    write(unit, '(A)') '      <Points>'
    write(unit, '(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    do i = 1, nvt
       write(unit, '(3(1X,ES24.16))') mesh%dcorvg(1, i), mesh%dcorvg(2, i), mesh%dcorvg(3, i)
    end do
    write(unit, '(A)') '        </DataArray>'
    write(unit, '(A)') '      </Points>'

    write(unit, '(A)') '      <Cells>'
    write(unit, '(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
    do i = 1, nel
       write(unit, '(*(1X,I10))') (int(mesh%kvert(j, i) - 1, kind=4), j = 1, nve)
    end do
    write(unit, '(A)') '        </DataArray>'

    write(unit, '(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
    offset = 0
    do i = 1, nel
       offset = offset + nve
       write(unit, '(1X,I10)') offset
    end do
    write(unit, '(A)') '        </DataArray>'

    write(unit, '(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
    do i = 1, nel
       write(unit, '(1X,I4)') 12
    end do
    write(unit, '(A)') '        </DataArray>'
    write(unit, '(A)') '      </Cells>'

    write(unit, '(A)') '    </Piece>'
    write(unit, '(A)') '  </UnstructuredGrid>'
    write(unit, '(A)') '</VTKFile>'

    close(unit)
    write(*, '(A,A)') 'VTU mesh written to ', trim(filename)
  end subroutine write_mesh_vtu

  subroutine compute_characteristic_size(mesh, config, characteristic_size)
    type(hex_mesh_type), intent(in) :: mesh
    type(MeshConfig), intent(in) :: config
    real(rk), intent(out) :: characteristic_size
    real(rk) :: bounds_min(3), bounds_max(3)
    real(rk) :: lengths(3), radial_span
    character(len=:), allocatable :: mesh_type

    call compute_mesh_bounds(mesh, bounds_min, bounds_max)
    lengths = bounds_max - bounds_min
    mesh_type = lowercase(trim(config%mesh_type))
    select case (trim(mesh_type))
    case ('box')
       characteristic_size = max(real(get_tolerance_floor(), rk), minval(max(lengths, 0.0_rk)))
    case ('hollowcylinder')
       radial_span = 0.0_rk
       if (config%cylinder%has_barrel_diameter .and. config%cylinder%has_inner_diameter) then
          radial_span = 0.5_rk * real(config%cylinder%barrel_diameter - config%cylinder%inner_diameter, rk)
       end if
       characteristic_size = max(real(get_tolerance_floor(), rk), min(lengths(3), radial_span))
    case default
       characteristic_size = max(real(get_tolerance_floor(), rk), minval(max(lengths, 0.0_rk)))
    end select
    if (characteristic_size <= 0.0_rk) characteristic_size = 1.0_rk
  end subroutine compute_characteristic_size

  subroutine compute_vertex_weights(distances, scale_factor, vertex_weights)
    real(rk), intent(in) :: distances(:)
    real(rk), intent(in) :: scale_factor
    real(rk), intent(out) :: vertex_weights(:)
    integer :: idx
    real(rk) :: d_scaled, kernel_value

    do idx = 1, size(distances)
      d_scaled = scale_factor * distances(idx)
      call kernel_function(d_scaled, kernel_value)
      kernel_value = min(kernel_value, 25.0_rk)
      vertex_weights(idx) = kernel_value**2.3_rk
    end do
  end subroutine compute_vertex_weights

  subroutine kernel_function(distance_value, response)
    real(rk), intent(inout) :: distance_value
    real(rk), intent(out) :: response
    real(rk) :: temp

    if (distance_value < 0.0_rk) distance_value = 2.5_rk * abs(distance_value)
    if (distance_value < 3.0_rk) then
       temp = 2.0_rk + 1.0_rk * (3.0_rk - distance_value)
    else
       temp = 2.0_rk - 0.2_rk * (distance_value - 3.0_rk)
    end if
    response = max(temp, 0.8_rk)
  end subroutine kernel_function

  subroutine compute_box_limits(mesh, config, start, finish)
    type(hex_mesh_type), intent(in) :: mesh
    type(MeshConfig), intent(in) :: config
    real(rk), intent(out) :: start(3)
    real(rk), intent(out) :: finish(3)

    call compute_mesh_bounds(mesh, start, finish)
    if (config%box%has_geometry_start) start = real(config%box%geometry_start, rk)
    if (config%box%has_geometry_length) finish = start + real(config%box%geometry_length, rk)
  end subroutine compute_box_limits

  subroutine compute_cylinder_parameters(mesh, config, center, inner_radius, outer_radius, axial_min, axial_max, &
       has_inner, has_outer, has_axial_min, has_axial_max)
    type(hex_mesh_type), intent(in) :: mesh
    type(MeshConfig), intent(in) :: config
    real(rk), intent(out) :: center(2)
    real(rk), intent(out) :: inner_radius, outer_radius
    real(rk), intent(out) :: axial_min, axial_max
    logical, intent(out) :: has_inner, has_outer, has_axial_min, has_axial_max
    real(rk) :: bounds_min(3), bounds_max(3)

    call compute_mesh_bounds(mesh, bounds_min, bounds_max)
    center(1) = 0.5_rk * (bounds_min(1) + bounds_max(1))
    center(2) = 0.5_rk * (bounds_min(2) + bounds_max(2))

    has_outer = config%cylinder%has_barrel_diameter
    if (has_outer) then
       outer_radius = 0.5_rk * real(config%cylinder%barrel_diameter, rk)
    else
       outer_radius = 0.0_rk
    end if
    has_inner = config%cylinder%has_inner_diameter
    if (has_inner) then
       inner_radius = 0.5_rk * real(config%cylinder%inner_diameter, rk)
    else
       inner_radius = 0.0_rk
    end if
    has_axial_min = config%cylinder%has_axial_start
    if (has_axial_min) then
       axial_min = real(config%cylinder%axial_start, rk)
    else
       axial_min = bounds_min(3)
    end if
    has_axial_max = has_axial_min .and. config%cylinder%has_barrel_length
    if (has_axial_max) then
       axial_max = real(config%cylinder%axial_start + config%cylinder%barrel_length, rk)
    else
       axial_max = bounds_max(3)
    end if
  end subroutine compute_cylinder_parameters

  integer function count_constraints(mask)
    integer, intent(in) :: mask
    integer :: value, bits
    value = mask
    bits = 0
    do while (value /= 0)
       if (iand(value, 1) /= 0) bits = bits + 1
       value = ishft(value, -1)
    end do
    count_constraints = bits
  end function count_constraints

  subroutine enforce_parametrization_constraints(constraint_flags, idx, current, px, py, pz, box_start, box_finish, &
       cyl_center, inner_radius, outer_radius, axial_min, axial_max, has_inner, has_outer, has_axial_min, has_axial_max)
    integer, intent(in) :: constraint_flags(:)
    integer, intent(in) :: idx
    real(rk), intent(in) :: current(3)
    real(rk), intent(inout) :: px, py, pz
    real(rk), intent(in) :: box_start(3), box_finish(3)
    real(rk), intent(in) :: cyl_center(2), inner_radius, outer_radius
    real(rk), intent(in) :: axial_min, axial_max
    logical, intent(in) :: has_inner, has_outer, has_axial_min, has_axial_max
    integer :: mask, combo_count
    logical :: has_xmin, has_xmax, has_ymin, has_ymax, has_zmin, has_zmax, has_inner_flag, has_outer_flag

    if (idx > size(constraint_flags)) return
    mask = constraint_flags(idx)
    combo_count = count_constraints(mask)
    if (combo_count <= 0) return
    if (combo_count >= 3) then
       px = current(1)
       py = current(2)
       pz = current(3)
       return
    end if

    has_xmin = iand(mask, constraint_x_min) /= 0
    has_xmax = iand(mask, constraint_x_max) /= 0
    has_ymin = iand(mask, constraint_y_min) /= 0
    has_ymax = iand(mask, constraint_y_max) /= 0
    has_zmin = iand(mask, constraint_z_min) /= 0
    has_zmax = iand(mask, constraint_z_max) /= 0
    has_inner_flag = iand(mask, constraint_cyl_inner) /= 0
    has_outer_flag = iand(mask, constraint_cyl_outer) /= 0

    select case (combo_count)
    case (1)
      call apply_single_constraint()
    case (2)
       if (apply_dual_constraint()) return
       call apply_single_constraint()
    end select
  contains
    subroutine apply_single_constraint()
      if (has_xmin) px = box_start(1)
      if (has_xmax) px = box_finish(1)
      if (has_ymin) py = box_start(2)
      if (has_ymax) py = box_finish(2)
      call apply_z_plane()
      if (has_inner_flag .and. has_inner) call project_to_radius(inner_radius)
      if (has_outer_flag .and. has_outer) call project_to_radius(outer_radius)
    end subroutine apply_single_constraint

    logical function apply_dual_constraint()
      apply_dual_constraint = .false.
      if ((has_xmin .or. has_xmax) .and. (has_ymin .or. has_ymax)) then
         if (has_xmin) px = box_start(1)
         if (has_xmax) px = box_finish(1)
         if (has_ymin) py = box_start(2)
         if (has_ymax) py = box_finish(2)
         apply_dual_constraint = .true.
         return
      end if
      if ((has_xmin .or. has_xmax) .and. (has_zmin .or. has_zmax)) then
         if (has_xmin) px = box_start(1)
         if (has_xmax) px = box_finish(1)
         call apply_z_plane()
         apply_dual_constraint = .true.
         return
      end if
      if ((has_ymin .or. has_ymax) .and. (has_zmin .or. has_zmax)) then
         if (has_ymin) py = box_start(2)
         if (has_ymax) py = box_finish(2)
         call apply_z_plane()
         apply_dual_constraint = .true.
         return
      end if
      if ((has_zmin .or. has_zmax) .and. (has_inner_flag .or. has_outer_flag)) then
         call apply_z_plane()
         if (has_inner_flag .and. has_inner) call project_to_radius(inner_radius)
         if (has_outer_flag .and. has_outer) call project_to_radius(outer_radius)
         apply_dual_constraint = .true.
         return
      end if
      if ((has_xmin .or. has_xmax .or. has_ymin .or. has_ymax) .and. (has_inner_flag .or. has_outer_flag)) then
         ! Intersections other than axial planes are unsupported; fall back to freezing
         px = current(1)
         py = current(2)
         pz = current(3)
         apply_dual_constraint = .true.
         return
      end if
    end function apply_dual_constraint

    subroutine project_to_radius(target_radius)
      real(rk), intent(in) :: target_radius
      real(rk) :: dx, dy, rad, scale
      dx = px - cyl_center(1)
      dy = py - cyl_center(2)
      rad = sqrt(dx * dx + dy * dy)
      if (rad <= 0.0_rk) then
         px = cyl_center(1) + target_radius
         py = cyl_center(2)
      else
         scale = target_radius / rad
         px = cyl_center(1) + dx * scale
         py = cyl_center(2) + dy * scale
      end if
    end subroutine project_to_radius

    subroutine apply_z_plane()
      if (has_zmin) then
         if (has_axial_min) then
            pz = axial_min
         else
            pz = box_start(3)
         end if
      end if
      if (has_zmax) then
         if (has_axial_max) then
            pz = axial_max
         else
            pz = box_finish(3)
         end if
      end if
    end subroutine apply_z_plane
  end subroutine enforce_parametrization_constraints

  subroutine compute_hexahedron_volumes(mesh, volumes)
    type(hex_mesh_type), intent(in) :: mesh
    real(rk), allocatable, intent(out) :: volumes(:)
    integer :: nel, v, cell_idx, face_idx
    integer :: v1, v2, v3, v4
    real(rk) :: cell_coords(3, 8)
    real(rk) :: centroid(3)
    real(rk) :: total

    if (.not.allocated(mesh%kvert)) then
       allocate(volumes(0))
       return
    end if

    if (size(mesh%kvert, 1) /= 8) error stop "Hex volume computation expects 8-node elements."
    nel = mesh%nel
    allocate(volumes(nel))
    volumes = 0.0_rk
    if (nel <= 0) return

    do cell_idx = 1, nel
       centroid = 0.0_rk
       do v = 1, 8
          cell_coords(:, v) = mesh%dcorvg(:, mesh%kvert(v, cell_idx))
          centroid = centroid + cell_coords(:, v)
       end do
       centroid = centroid / 8.0_rk
       total = 0.0_rk
       do face_idx = 1, size(hexahedron_faces, 2)
          v1 = hexahedron_faces(1, face_idx)
          v2 = hexahedron_faces(2, face_idx)
          v3 = hexahedron_faces(3, face_idx)
          v4 = hexahedron_faces(4, face_idx)
          total = total + tetrahedron_volume(centroid, cell_coords(:, v1), cell_coords(:, v2), cell_coords(:, v3))
          total = total + tetrahedron_volume(centroid, cell_coords(:, v1), cell_coords(:, v3), cell_coords(:, v4))
       end do
       volumes(cell_idx) = total
    end do
  end subroutine compute_hexahedron_volumes

  real(rk) function tetrahedron_volume(a, b, c, d) result(volume)
    real(rk), intent(in) :: a(3), b(3), c(3), d(3)
    real(rk) :: ab(3), ac(3), ad(3), cross_prod(3)

    ab = b - a
    ac = c - a
    ad = d - a
    cross_prod(1) = ac(2) * ad(3) - ac(3) * ad(2)
    cross_prod(2) = ac(3) * ad(1) - ac(1) * ad(3)
    cross_prod(3) = ac(1) * ad(2) - ac(2) * ad(1)
    volume = abs(ab(1) * cross_prod(1) + ab(2) * cross_prod(2) + ab(3) * cross_prod(3)) / 6.0_rk
  end function tetrahedron_volume

  subroutine build_edge_to_element(mesh)
    type(hex_mesh_type), intent(inout) :: mesh
    integer :: nel, elem_idx, edge_idx
    integer :: v1, v2, edge_id

    if (.not.allocated(mesh%kvert)) return
    nel = mesh%nel
    if (nel <= 0) return

    if (allocated(mesh%kedge)) deallocate(mesh%kedge)
    allocate(mesh%kedge(12, nel))
    mesh%kedge = 0

    mesh%nedge = 0
    mesh%edge_capacity = max(1, 6 * nel)
    if (allocated(mesh%edge_nodes)) deallocate(mesh%edge_nodes)
    allocate(mesh%edge_nodes(2, mesh%edge_capacity))

    do elem_idx = 1, nel
       do edge_idx = 1, size(hexahedron_edges, 2)
          v1 = int(mesh%kvert(hexahedron_edges(1, edge_idx), elem_idx))
          v2 = int(mesh%kvert(hexahedron_edges(2, edge_idx), elem_idx))
          if (v1 > v2) then
             call swap_int(v1, v2)
          end if
          edge_id = find_edge_id(mesh, v1, v2)
          mesh%kedge(edge_idx, elem_idx) = edge_id
       end do
    end do

    call trim_edge_storage(mesh)
  end subroutine build_edge_to_element

  subroutine broadcast_mesh(mesh, mpi_rank, mpi_size)
    type(hex_mesh_type), intent(inout) :: mesh
    integer, intent(in) :: mpi_rank, mpi_size
    integer :: header(6)
    integer :: mpi_err

    if (mpi_rank == 0) then
       header = [mesh%nel, mesh%nvt, mesh%nbct, mesh%nve, mesh%nee, mesh%nae]
    else
       header = 0
    end if
    call MPI_Bcast(header, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    if (mpi_rank /= 0) then
       mesh%nel = header(1)
       mesh%nvt = header(2)
       mesh%nbct = header(3)
       mesh%nve = header(4)
       mesh%nee = header(5)
       mesh%nae = header(6)
       if (.not.allocated(mesh%dcorvg)) allocate(mesh%dcorvg(3, mesh%nvt))
       if (.not.allocated(mesh%kvert)) allocate(mesh%kvert(mesh%nve, mesh%nel))
       if (.not.allocated(mesh%knpr)) allocate(mesh%knpr(mesh%nvt))
    end if
    if (mesh%nvt > 0) then
       call MPI_Bcast(mesh%dcorvg, 3 * mesh%nvt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)
       call MPI_Bcast(mesh%knpr, mesh%nvt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    end if
    if (mesh%nel > 0) then
       call MPI_Bcast(mesh%kvert, mesh%nve * mesh%nel, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    end if
  end subroutine broadcast_mesh

  subroutine broadcast_mesh_coordinates(mesh, mpi_rank)
    type(hex_mesh_type), intent(inout) :: mesh
    integer, intent(in) :: mpi_rank
    integer :: mpi_err

    if (.not.allocated(mesh%dcorvg)) return
    if (mesh%nvt <= 0) return
    call MPI_Bcast(mesh%dcorvg, 3 * mesh%nvt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)
  end subroutine broadcast_mesh_coordinates

  integer function find_edge_id(mesh, v1, v2) result(edge_id)
    type(hex_mesh_type), intent(inout) :: mesh
    integer, intent(in) :: v1, v2
    integer :: idx

    do idx = 1, mesh%nedge
       if (mesh%edge_nodes(1, idx) == v1 .and. mesh%edge_nodes(2, idx) == v2) then
          edge_id = idx
          return
       end if
    end do

    mesh%nedge = mesh%nedge + 1
    if (mesh%nedge > mesh%edge_capacity) then
       call grow_edge_storage(mesh, max(2 * mesh%edge_capacity, mesh%edge_capacity + 1))
    end if
    mesh%edge_nodes(1, mesh%nedge) = v1
    mesh%edge_nodes(2, mesh%nedge) = v2
    edge_id = mesh%nedge
  end function find_edge_id

  subroutine grow_edge_storage(mesh, new_capacity)
    type(hex_mesh_type), intent(inout) :: mesh
    integer, intent(in) :: new_capacity
    integer, allocatable :: temp(:,:)

    allocate(temp(2, new_capacity))
    temp = 0
    if (mesh%nedge > 0) temp(:, 1:mesh%nedge) = mesh%edge_nodes(:, 1:mesh%nedge)
    if (allocated(mesh%edge_nodes)) deallocate(mesh%edge_nodes)
    call move_alloc(temp, mesh%edge_nodes)
    mesh%edge_capacity = new_capacity
  end subroutine grow_edge_storage

  subroutine trim_edge_storage(mesh)
    type(hex_mesh_type), intent(inout) :: mesh
    integer, allocatable :: temp(:,:)

    if (mesh%nedge <= 0) then
       mesh%edge_capacity = 0
       if (allocated(mesh%edge_nodes)) deallocate(mesh%edge_nodes)
       return
    end if

    if (mesh%nedge == mesh%edge_capacity) return
    allocate(temp(2, mesh%nedge))
    temp = mesh%edge_nodes(:, 1:mesh%nedge)
    if (allocated(mesh%edge_nodes)) deallocate(mesh%edge_nodes)
    call move_alloc(temp, mesh%edge_nodes)
    mesh%edge_capacity = mesh%nedge
  end subroutine trim_edge_storage

  subroutine swap_int(a, b)
    integer, intent(inout) :: a, b
    integer :: tmp
    tmp = a
    a = b
    b = tmp
  end subroutine swap_int

  real(rk) function compute_min_edge_length(mesh) result(min_len)
    type(hex_mesh_type), intent(in) :: mesh
    integer :: nel, edge_idx, elem_idx
    integer :: vid_a, vid_b
    real(rk) :: diff(3), edge_len

    min_len = 0.0_rk
    if (.not.allocated(mesh%dcorvg)) return
    if (.not.allocated(mesh%kvert)) return
    nel = size(mesh%kvert, 2)
    if (nel <= 0) return

    min_len = huge(1.0_rk)
    do elem_idx = 1, nel
       do edge_idx = 1, size(hexahedron_edges, 2)
          vid_a = int(mesh%kvert(hexahedron_edges(1, edge_idx), elem_idx))
          vid_b = int(mesh%kvert(hexahedron_edges(2, edge_idx), elem_idx))
          if (vid_a <= 0 .or. vid_a > size(mesh%dcorvg, 2)) cycle
          if (vid_b <= 0 .or. vid_b > size(mesh%dcorvg, 2)) cycle
          diff = mesh%dcorvg(:, vid_a) - mesh%dcorvg(:, vid_b)
          edge_len = sqrt(sum(diff * diff))
          if (edge_len > 0.0_rk .and. edge_len < min_len) min_len = edge_len
       end do
    end do

    if (min_len == huge(1.0_rk)) min_len = 0.0_rk
  end function compute_min_edge_length

  subroutine identify_boundary_nodes(mesh, boundary_mask)
    type(hex_mesh_type), intent(in) :: mesh
    logical, allocatable, intent(out) :: boundary_mask(:)
    integer :: nel, nvt, total_faces
    integer :: face_idx, elem_idx, face_id, corner_idx
    integer, allocatable :: face_vertices(:, :)
    integer, allocatable :: face_keys(:, :)
    integer, allocatable :: face_order(:)
    logical, allocatable :: boundary_face(:)
    integer :: group_start, group_end
    integer :: idx

    nvt = mesh%nvt
    if (nvt <= 0) then
       allocate(boundary_mask(0))
       return
    end if

    allocate(boundary_mask(nvt))
    boundary_mask = .false.
    if (.not.allocated(mesh%kvert)) return
    nel = size(mesh%kvert, 2)
    if (nel <= 0) return

    total_faces = 6 * nel
    allocate(face_vertices(4, total_faces))
    allocate(face_keys(4, total_faces))
    face_idx = 0
    do elem_idx = 1, nel
       do face_id = 1, 6
          face_idx = face_idx + 1
          do corner_idx = 1, 4
             face_vertices(corner_idx, face_idx) = int(mesh%kvert(hexahedron_faces(corner_idx, face_id), elem_idx))
             face_keys(corner_idx, face_idx) = face_vertices(corner_idx, face_idx)
          end do
          call sort_face_vertices(face_keys(:, face_idx))
       end do
    end do

    allocate(face_order(total_faces))
    do idx = 1, total_faces
       face_order(idx) = idx
    end do
    call sort_face_order(face_keys, face_order)

    allocate(boundary_face(total_faces))
    boundary_face = .false.
    face_idx = 1
    do while (face_idx <= total_faces)
       group_start = face_idx
       group_end = face_idx
       do while (group_end < total_faces)
          if (.not.faces_equal(face_keys, face_order(group_end), face_order(group_end + 1))) exit
          group_end = group_end + 1
       end do
       if (group_start == group_end) boundary_face(face_order(group_start)) = .true.
       face_idx = group_end + 1
    end do

    do face_idx = 1, total_faces
       if (.not.boundary_face(face_idx)) cycle
       do corner_idx = 1, 4
          if (face_vertices(corner_idx, face_idx) >= 1 .and. face_vertices(corner_idx, face_idx) <= nvt) then
             boundary_mask(face_vertices(corner_idx, face_idx)) = .true.
          end if
       end do
    end do

    deallocate(face_vertices, face_keys, face_order, boundary_face)
  end subroutine identify_boundary_nodes

  subroutine sort_face_vertices(face)
    integer, intent(inout) :: face(4)
    integer :: i, j, temp

    do i = 2, 4
       j = i
       do while (j > 1 .and. face(j-1) > face(j))
          temp = face(j-1)
          face(j-1) = face(j)
          face(j) = temp
          j = j - 1
       end do
    end do
  end subroutine sort_face_vertices

  subroutine sort_face_order(face_keys, face_order)
    integer, intent(in) :: face_keys(:, :)
    integer, intent(inout) :: face_order(:)
    integer :: i, j, key

    do i = 2, size(face_order)
       key = face_order(i)
       j = i - 1
       do while (j >= 1 .and. compare_faces(face_keys, face_order(j), key) > 0)
          face_order(j + 1) = face_order(j)
          j = j - 1
       end do
       face_order(j + 1) = key
    end do
  end subroutine sort_face_order

  integer function compare_faces(face_keys, idx_a, idx_b)
    integer, intent(in) :: face_keys(:, :)
    integer, intent(in) :: idx_a, idx_b
    integer :: k

    compare_faces = 0
    do k = 1, min(size(face_keys, 1), 4)
       if (face_keys(k, idx_a) < face_keys(k, idx_b)) then
          compare_faces = -1
          return
       else if (face_keys(k, idx_a) > face_keys(k, idx_b)) then
          compare_faces = 1
          return
       end if
    end do
  end function compare_faces

  logical function faces_equal(face_keys, idx_a, idx_b) result(equal)
    integer, intent(in) :: face_keys(:, :)
    integer, intent(in) :: idx_a, idx_b
    integer :: k

    equal = .true.
    do k = 1, min(size(face_keys, 1), 4)
       if (face_keys(k, idx_a) /= face_keys(k, idx_b)) then
          equal = .false.
          return
       end if
    end do
  end function faces_equal

  subroutine report_boundary_assignment_status(constraint_flags, boundary_mask)
    integer, intent(in) :: constraint_flags(:)
    logical, intent(in) :: boundary_mask(:)
    integer :: boundary_total, missing

    if (size(boundary_mask) == 0) return
    boundary_total = count(boundary_mask)
    if (boundary_total <= 0) return
    missing = count(boundary_mask .and. (constraint_flags == constraint_none))
    if (missing > 0) then
       write(*, '(A,I0)') 'Boundary vertices lacking parametrization: ', missing
    else
       write(*, '(A)') 'All boundary vertices assigned to parametrizations.'
    end if
  end subroutine report_boundary_assignment_status

  subroutine report_constraint_combinations(constraint_flags)
    integer, intent(in) :: constraint_flags(:)
    integer, allocatable :: combos(:)
    integer, allocatable :: counts(:)
    integer :: combo_count, idx
    character(len=128) :: label

    call collect_combinations(constraint_flags, combos, counts, combo_count)
    if (combo_count <= 0) then
       write(*, '(A)') 'No parametrized vertices detected.'
       return
    end if

    write(*, '(A)') 'Combination assignments:'
    do idx = 1, combo_count
       if (counts(idx) <= 0) cycle
       call format_constraint_combo(combos(idx), label)
       write(*, '(A,1X,I0)') '  [' // trim(label) // ']:', counts(idx)
    end do

    if (allocated(combos)) deallocate(combos)
    if (allocated(counts)) deallocate(counts)
  end subroutine report_constraint_combinations

  subroutine collect_combinations(constraint_flags, combos, counts, combo_count)
    integer, intent(in) :: constraint_flags(:)
    integer, allocatable, intent(out) :: combos(:)
    integer, allocatable, intent(out) :: counts(:)
    integer, intent(out) :: combo_count
    integer :: idx, pos

    combo_count = 0
    allocate(combos(0))
    allocate(counts(0))
    do idx = 1, size(constraint_flags)
       if (constraint_flags(idx) == 0) cycle
       pos = find_combo(constraint_flags(idx), combos, combo_count)
       if (pos > 0) then
          counts(pos) = counts(pos) + 1
       else
          call append_combo(constraint_flags(idx), combos, counts, combo_count)
       end if
    end do
    call sort_combos(combos, counts, combo_count)
  end subroutine collect_combinations

  integer function find_combo(value, combos, combo_count) result(pos)
    integer, intent(in) :: value
    integer, intent(in) :: combos(:)
    integer, intent(in) :: combo_count
    integer :: idx

    pos = 0
    do idx = 1, combo_count
       if (combos(idx) == value) then
          pos = idx
          return
       end if
    end do
  end function find_combo

  subroutine append_combo(value, combos, counts, combo_count)
    integer, intent(in) :: value
    integer, allocatable, intent(inout) :: combos(:)
    integer, allocatable, intent(inout) :: counts(:)
    integer, intent(inout) :: combo_count
    integer, allocatable :: new_combos(:), new_counts(:)

    combo_count = combo_count + 1
    allocate(new_combos(combo_count))
    allocate(new_counts(combo_count))
    if (combo_count > 1) then
       new_combos(1:combo_count-1) = combos
       new_counts(1:combo_count-1) = counts
    end if
    new_combos(combo_count) = value
    new_counts(combo_count) = 1
    if (allocated(combos)) deallocate(combos)
    if (allocated(counts)) deallocate(counts)
    call move_alloc(new_combos, combos)
    call move_alloc(new_counts, counts)
  end subroutine append_combo

  subroutine sort_combos(combos, counts, combo_count)
    integer, intent(inout) :: combos(:)
    integer, intent(inout) :: counts(:)
    integer, intent(in) :: combo_count
    integer :: i, j, tmp_val, tmp_count

    do i = 1, combo_count - 1
       do j = i + 1, combo_count
          if (counts(j) > counts(i)) then
             tmp_val = combos(i)
             combos(i) = combos(j)
             combos(j) = tmp_val
             tmp_count = counts(i)
             counts(i) = counts(j)
             counts(j) = tmp_count
          end if
       end do
    end do
  end subroutine sort_combos

  subroutine format_constraint_combo(value, label)
    integer, intent(in) :: value
    character(len=*), intent(out) :: label
    logical :: first
    integer :: idx

    label = ''
    first = .true.
    do idx = 1, size(constraint_masks)
       if (iand(value, constraint_masks(idx)) /= 0) then
          if (.not.first) label = trim(label) // '|'
          label = trim(label) // trim(constraint_names(idx))
          first = .false.
       end if
    end do
    if (first) label = 'none'
  end subroutine format_constraint_combo

  subroutine collect_nodes_from_mask(mask, list)
    logical, intent(in) :: mask(:)
    integer, allocatable, intent(out) :: list(:)
    integer :: count_value, idx, pos

    count_value = count(mask)
    if (allocated(list)) deallocate(list)
    allocate(list(count_value))
    if (count_value == 0) return
    pos = 0
    do idx = 1, size(mask)
       if (.not.mask(idx)) cycle
       pos = pos + 1
       list(pos) = idx
    end do
  end subroutine collect_nodes_from_mask

  subroutine build_constraint_lists(constraint_flags, lists, boundary_mask)
    integer, intent(in) :: constraint_flags(:)
    type(vertex_constraint_lists), intent(out) :: lists
    logical, intent(in), optional :: boundary_mask(:)

    call fill_constraint_list(constraint_flags, constraint_x_min, lists%x_min)
    call fill_constraint_list(constraint_flags, constraint_x_max, lists%x_max)
    call fill_constraint_list(constraint_flags, constraint_y_min, lists%y_min)
    call fill_constraint_list(constraint_flags, constraint_y_max, lists%y_max)
    call fill_constraint_list(constraint_flags, constraint_z_min, lists%z_min)
    call fill_constraint_list(constraint_flags, constraint_z_max, lists%z_max)
    call fill_constraint_list(constraint_flags, constraint_cyl_inner, lists%cyl_inner)
    call fill_constraint_list(constraint_flags, constraint_cyl_outer, lists%cyl_outer)
    if (present(boundary_mask)) then
       call collect_nodes_from_mask(boundary_mask, lists%boundary)
    else
       if (allocated(lists%boundary)) deallocate(lists%boundary)
       allocate(lists%boundary(0))
    end if
  end subroutine build_constraint_lists

  subroutine fill_constraint_list(constraint_flags, bitmask, list)
    integer, intent(in) :: constraint_flags(:)
    integer, intent(in) :: bitmask
    integer, allocatable, intent(out) :: list(:)
    integer :: count_value, idx, pos

    if (allocated(list)) deallocate(list)
    count_value = count(iand(constraint_flags, bitmask) /= 0)
    allocate(list(count_value))
    if (count_value == 0) return
    pos = 0
    do idx = 1, size(constraint_flags)
       if (iand(constraint_flags(idx), bitmask) == 0) cycle
       pos = pos + 1
       list(pos) = idx
    end do
  end subroutine fill_constraint_list

  subroutine release_constraint_lists(lists)
    type(vertex_constraint_lists), intent(inout) :: lists

    if (allocated(lists%x_min)) deallocate(lists%x_min)
    if (allocated(lists%x_max)) deallocate(lists%x_max)
    if (allocated(lists%y_min)) deallocate(lists%y_min)
    if (allocated(lists%y_max)) deallocate(lists%y_max)
    if (allocated(lists%z_min)) deallocate(lists%z_min)
    if (allocated(lists%z_max)) deallocate(lists%z_max)
    if (allocated(lists%cyl_inner)) deallocate(lists%cyl_inner)
    if (allocated(lists%cyl_outer)) deallocate(lists%cyl_outer)
    if (allocated(lists%boundary)) deallocate(lists%boundary)
  end subroutine release_constraint_lists

end module def_mod
