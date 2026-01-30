module bc_treatment
    use iso_c_binding, only: c_double, c_int
    implicit none
    private

    integer, parameter :: hexahedron_edges(2, 12) = reshape([ &
        1, 2, 2, 3, 3, 4, 4, 1, 5, 6, 6, 7, 7, 8, 8, 5, 1, 5, 2, 6, 3, 7, 4, 8], [2, 12])
    integer, parameter :: hexahedron_faces(4, 6) = reshape([ &
        1, 2, 3, 4, &
        5, 6, 7, 8, &
        1, 2, 6, 5, &
        2, 3, 7, 6, &
        3, 4, 8, 7, &
        4, 1, 5, 8], [4, 6])

    type, public :: FaceList
        integer :: count = 0
        integer, allocatable :: elements(:)
        integer, allocatable :: face_ids(:)
    end type FaceList

    type, public :: BoxMeshConfig
        real(c_double) :: geometry_start(3) = 0.0_c_double
        real(c_double) :: geometry_length(3) = 0.0_c_double
        real(c_double) :: spacing(3) = 0.0_c_double
        logical :: has_geometry_start = .false.
        logical :: has_geometry_length = .false.
        logical :: has_spacing = .false.
    end type BoxMeshConfig

    type, public :: HollowCylinderMeshConfig
        real(c_double) :: barrel_diameter = 0.0_c_double
        real(c_double) :: inner_diameter = 0.0_c_double
        real(c_double) :: barrel_length = 0.0_c_double
        real(c_double) :: axial_start = 0.0_c_double
        real(c_double) :: spacing(3) = 0.0_c_double
        logical :: has_barrel_diameter = .false.
        logical :: has_inner_diameter = .false.
        logical :: has_barrel_length = .false.
        logical :: has_axial_start = .false.
        logical :: has_spacing = .false.
    end type HollowCylinderMeshConfig

    type, public :: MeshConfig
        character(len=64) :: mesh_type = "UNSPECIFIED"
        logical :: loaded = .false.
        character(len=512) :: source_file = ""
        type(BoxMeshConfig) :: box
        type(HollowCylinderMeshConfig) :: cylinder
    end type MeshConfig

    type, public :: HollowCylinderBoundaryClassification
        real(c_double) :: min_edge_length = 0.0_c_double
        real(c_double) :: tolerance = 0.0_c_double
        integer :: total_boundary_nodes = 0
        integer, allocatable :: cyl_outer_nodes(:)
        integer, allocatable :: cyl_inner_nodes(:)
        integer, allocatable :: axial_min_nodes(:)
        integer, allocatable :: axial_max_nodes(:)
        integer, allocatable :: inner_wall_nodes(:)
        type(FaceList) :: all_boundary_faces
        type(FaceList) :: cyl_outer_faces
        type(FaceList) :: cyl_inner_faces
        type(FaceList) :: axial_min_faces
        type(FaceList) :: axial_max_faces
    end type HollowCylinderBoundaryClassification

    type, public :: BoxBoundaryClassification
        real(c_double) :: min_edge_length = 0.0_c_double
        real(c_double) :: tolerance = 0.0_c_double
        integer :: total_boundary_nodes = 0
        integer, allocatable :: x_min_nodes(:)
        integer, allocatable :: x_max_nodes(:)
        integer, allocatable :: y_min_nodes(:)
        integer, allocatable :: y_max_nodes(:)
        integer, allocatable :: z_min_nodes(:)
        integer, allocatable :: z_max_nodes(:)
        integer, allocatable :: inner_wall_nodes(:)
        type(FaceList) :: all_boundary_faces
        type(FaceList) :: x_min_faces
        type(FaceList) :: x_max_faces
        type(FaceList) :: y_min_faces
        type(FaceList) :: y_max_faces
        type(FaceList) :: z_min_faces
        type(FaceList) :: z_max_faces
    end type BoxBoundaryClassification

    public :: initialize_mesh_config, load_mesh_config, log_mesh_config
    public :: recompute_knpr_from_connectivity
    public :: classify_hollowcylinder_boundaries, classify_box_boundaries

contains
    subroutine initialize_mesh_config(config)
        type(MeshConfig), intent(out) :: config

        config%mesh_type = "UNSPECIFIED"
        config%loaded = .false.
        config%source_file = ""
        config%box = BoxMeshConfig()
        config%cylinder = HollowCylinderMeshConfig()
    end subroutine initialize_mesh_config

    subroutine load_mesh_config(filename, config)
        character(len=*), intent(in) :: filename
        type(MeshConfig), intent(inout) :: config

        logical :: exists
        integer :: unit, ios, eq_pos, line_len
        character(len=512) :: raw_line, line, section, key, value

        inquire(file=filename, exist=exists)
        if (.not. exists) then
            write(*, '(A)') "Mesh configuration file not found: "//trim(filename)
            return
        end if

        call initialize_mesh_config(config)
        config%source_file = trim(filename)
        section = ""

        open(newunit=unit, file=filename, status="old", action="read")
        do
            read(unit, '(A)', iostat=ios) raw_line
            if (ios /= 0) exit
            line = adjustl(trim(raw_line))
            line_len = len_trim(line)
            if (line_len == 0) cycle
            if (line(1:1) == "!" .or. line(1:1) == "#") cycle
            if (line(1:1) == "[") then
                call parse_section_name(line, section)
                cycle
            end if
            eq_pos = index(line, "=")
            if (eq_pos <= 1 .or. eq_pos >= line_len) cycle
            key = trim(adjustl(line(:eq_pos - 1)))
            value = trim(adjustl(line(eq_pos + 1:)))
            call assign_config_value(section, key, value, config)
        end do
        close(unit)
        config%loaded = .true.
    end subroutine load_mesh_config

    subroutine log_mesh_config(config)
        type(MeshConfig), intent(in) :: config
        character(len=32) :: type_token

        if (.not. config%loaded) then
            return
        end if

        write(*,'(A)') "Mesh configuration summary:"
        write(*,'(A,1X,A)') "  Source:", trim(config%source_file)
        write(*,'(A,1X,A)') "  HexMesher:", trim(config%mesh_type)
        type_token = lowercase(trim(config%mesh_type))
        select case (trim(type_token))
        case ("box")
            if (config%box%has_geometry_start) &
                write(*,'(A,3(1X,ES12.5))') "  geometryStart:", config%box%geometry_start
            if (config%box%has_geometry_length) &
                write(*,'(A,3(1X,ES12.5))') "  geometryLength:", config%box%geometry_length
            if (config%box%has_spacing) &
                write(*,'(A,3(1X,ES12.5))') "  sEl(x,y,z):", config%box%spacing
        case ("hollowcylinder")
            if (config%cylinder%has_barrel_diameter) &
                write(*,'(A,1X,ES12.5)') "  BarrelDiameter:", config%cylinder%barrel_diameter
            if (config%cylinder%has_inner_diameter) &
                write(*,'(A,1X,ES12.5)') "  InnerDiameter :", config%cylinder%inner_diameter
            if (config%cylinder%has_barrel_length) &
                write(*,'(A,1X,ES12.5)') "  BarrelLength  :", config%cylinder%barrel_length
            if (config%cylinder%has_axial_start) &
                write(*,'(A,1X,ES12.5)') "  AxialStartPos :", config%cylinder%axial_start
            if (config%cylinder%has_spacing) &
                write(*,'(A,3(1X,ES12.5))') "  sEl(tan,rad,ax):", config%cylinder%spacing
        end select
    end subroutine log_mesh_config

    subroutine assign_config_value(section, key, value, config)
        character(len=*), intent(in) :: section, key, value
        type(MeshConfig), intent(inout) :: config
        character(len=128) :: section_id, key_id
        real(c_double) :: temp_vec(3)
        logical :: parsed

        section_id = lowercase(trim(section))
        key_id = lowercase(trim(key))

        select case (section_id)
        case ("e3dsimulationsettings")
            select case (key_id)
            case ("hexmesher")
                call set_mesh_type(config, value)
            case ("sel_x")
                parsed = try_parse_real(value, config%box%spacing(1))
                if (parsed) config%box%has_spacing = .true.
            case ("sel_y")
                parsed = try_parse_real(value, config%box%spacing(2))
                if (parsed) config%box%has_spacing = .true.
            case ("sel_z")
                parsed = try_parse_real(value, config%box%spacing(3))
                if (parsed) config%box%has_spacing = .true.
            case ("sel_tangential")
                parsed = try_parse_real(value, config%cylinder%spacing(1))
                if (parsed) config%cylinder%has_spacing = .true.
            case ("sel_radial")
                parsed = try_parse_real(value, config%cylinder%spacing(2))
                if (parsed) config%cylinder%has_spacing = .true.
            case ("sel_axial")
                parsed = try_parse_real(value, config%cylinder%spacing(3))
                if (parsed) config%cylinder%has_spacing = .true.
            end select
        case ("e3dgeometrydata/machine")
            select case (key_id)
            case ("geometrystart")
                parsed = try_parse_real_vector(value, temp_vec)
                if (parsed) then
                    config%box%geometry_start = temp_vec
                    config%box%has_geometry_start = .true.
                end if
            case ("geometrylength")
                parsed = try_parse_real_vector(value, temp_vec)
                if (parsed) then
                    config%box%geometry_length = temp_vec
                    config%box%has_geometry_length = .true.
                end if
            case ("barreldiameter")
                parsed = try_parse_real(value, config%cylinder%barrel_diameter)
                if (parsed) config%cylinder%has_barrel_diameter = .true.
            case ("innerdiameter")
                parsed = try_parse_real(value, config%cylinder%inner_diameter)
                if (parsed) config%cylinder%has_inner_diameter = .true.
            case ("barrellength")
                parsed = try_parse_real(value, config%cylinder%barrel_length)
                if (parsed) config%cylinder%has_barrel_length = .true.
            case ("axialstartposition")
                parsed = try_parse_real(value, config%cylinder%axial_start)
                if (parsed) config%cylinder%has_axial_start = .true.
            end select
        end select
    end subroutine assign_config_value

    subroutine set_mesh_type(config, raw_value)
        type(MeshConfig), intent(inout) :: config
        character(len=*), intent(in) :: raw_value
        character(len=64) :: token

        token = lowercase(trim(raw_value))
        select case (trim(token))
        case ("box")
            config%mesh_type = "Box"
        case ("hollowcylinder")
            config%mesh_type = "HollowCylinder"
        case default
            config%mesh_type = trim(raw_value)
        end select
    end subroutine set_mesh_type

    logical function try_parse_real(str, result_value)
        character(len=*), intent(in) :: str
        real(c_double), intent(out) :: result_value
        integer :: ios

        read(str, *, iostat=ios) result_value
        try_parse_real = (ios == 0)
    end function try_parse_real

    logical function try_parse_real_vector(str, values)
        character(len=*), intent(in) :: str
        real(c_double), intent(out) :: values(:)
        character(len=512) :: sanitized
        integer :: ios, idx

        sanitized = str
        do idx = 1, len(sanitized)
            if (sanitized(idx:idx) == ",") sanitized(idx:idx) = " "
        end do
        read(sanitized, *, iostat=ios) values
        try_parse_real_vector = (ios == 0)
    end function try_parse_real_vector

    subroutine parse_section_name(line, section)
        character(len=*), intent(in) :: line
        character(len=*), intent(inout) :: section
        integer :: close_idx

        close_idx = index(line, "]")
        if (close_idx > 2) then
            section = trim(line(2:close_idx - 1))
        else
            section = ""
        end if
    end subroutine parse_section_name

    pure function lowercase(str) result(out)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: out
        integer :: idx, code

        do idx = 1, len(str)
            code = iachar(str(idx:idx))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                out(idx:idx) = achar(code + 32)
            else
                out(idx:idx) = str(idx:idx)
            end if
        end do
    end function lowercase

    subroutine recompute_knpr_from_connectivity(kvert, knpr, boundary_faces)
        integer(c_int), intent(in) :: kvert(:, :)
        integer(c_int), intent(inout) :: knpr(:)
        type(FaceList), intent(out), optional :: boundary_faces

        integer :: nel, nve, total_faces
        integer :: face_idx, elem_idx, face_id, corner_idx
        integer :: vid, group_start, group_end
        integer, allocatable :: face_vertices(:, :)
        integer, allocatable :: face_keys(:, :)
        integer, allocatable :: face_order(:)
        logical, allocatable :: boundary_face(:)
        integer, allocatable :: face_elem(:), face_local(:)
        integer :: boundary_count, pos

        nve = size(kvert, 1)
        if (nve /= 8) then
            error stop "KNPR boundary tagging currently assumes 8-node hexahedra."
        end if

        nel = size(kvert, 2)
        if (nel == 0 .or. size(knpr) == 0) then
            if (size(knpr) > 0) knpr = 0
            return
        end if

        total_faces = 6 * nel
        knpr = 0

        allocate(face_vertices(4, total_faces))
        allocate(face_keys(4, total_faces))
        allocate(face_elem(total_faces))
        allocate(face_local(total_faces))
        face_idx = 0
        do elem_idx = 1, nel
            do face_id = 1, 6
                face_idx = face_idx + 1
                do corner_idx = 1, 4
                    vid = int(kvert(hexahedron_faces(corner_idx, face_id), elem_idx))
                    face_vertices(corner_idx, face_idx) = vid
                    face_keys(corner_idx, face_idx) = vid
                end do
                face_elem(face_idx) = elem_idx
                face_local(face_idx) = face_id
                call sort_face_vertices(face_keys(:, face_idx))
            end do
        end do

        allocate(face_order(total_faces))
        do face_idx = 1, total_faces
            face_order(face_idx) = face_idx
        end do
        call sort_faces_by_keys(face_keys, face_order)

        allocate(boundary_face(total_faces))
        boundary_face = .false.
        face_idx = 1
        do while (face_idx <= total_faces)
            group_start = face_idx
            group_end = face_idx
            do while (group_end < total_faces)
                if (.not. faces_equal(face_keys, face_order(group_end), face_order(group_end + 1))) exit
                group_end = group_end + 1
            end do
            if (group_end == group_start) boundary_face(face_order(group_start)) = .true.
            face_idx = group_end + 1
        end do

        do face_idx = 1, total_faces
            if (.not. boundary_face(face_idx)) cycle
            do corner_idx = 1, 4
                vid = face_vertices(corner_idx, face_idx)
                knpr(vid) = 1_c_int
            end do
        end do

        if (present(boundary_faces)) then
            boundary_count = count(boundary_face)
            call allocate_face_list(boundary_faces, boundary_count)
            if (boundary_count > 0) then
                pos = 0
                do face_idx = 1, total_faces
                    if (.not. boundary_face(face_idx)) cycle
                    pos = pos + 1
                    boundary_faces%elements(pos) = face_elem(face_idx)
                    boundary_faces%face_ids(pos) = face_local(face_idx)
                end do
            end if
        end if

        deallocate(face_vertices, face_keys, face_order, boundary_face, face_elem, face_local)
    end subroutine recompute_knpr_from_connectivity

    subroutine classify_hollowcylinder_boundaries(config, coords, kvert, knpr, boundary_faces, classification)
        type(MeshConfig), intent(in) :: config
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: kvert(:, :)
        integer(c_int), intent(in) :: knpr(:)
        type(FaceList), intent(in) :: boundary_faces
        type(HollowCylinderBoundaryClassification), intent(out) :: classification

        real(c_double) :: outer_radius, inner_radius, axial_min, axial_max
        logical :: have_outer, have_inner, have_axial_min, have_axial_max
        integer :: face_count
        integer :: face_idx, node_idx
        integer :: node_ids(4)
        real(c_double) :: radial, axial
        logical, allocatable :: face_is_outer(:), face_is_inner(:)
        logical, allocatable :: face_is_axial_min(:), face_is_axial_max(:)

        classification%min_edge_length = compute_min_edge_length(coords, kvert)
        classification%tolerance = 0.90_c_double * classification%min_edge_length
        classification%total_boundary_nodes = count(knpr == 1_c_int)
        call allocate_zero_length_list(classification%cyl_outer_nodes)
        call allocate_zero_length_list(classification%cyl_inner_nodes)
        call allocate_zero_length_list(classification%axial_min_nodes)
        call allocate_zero_length_list(classification%axial_max_nodes)
        call allocate_zero_length_list(classification%inner_wall_nodes)
        call init_face_list(classification%cyl_outer_faces)
        call init_face_list(classification%cyl_inner_faces)
        call init_face_list(classification%axial_min_faces)
        call init_face_list(classification%axial_max_faces)
        call init_face_list(classification%all_boundary_faces)

        if (classification%total_boundary_nodes == 0) return

        have_outer = config%cylinder%has_barrel_diameter
        if (have_outer) outer_radius = 0.5_c_double * config%cylinder%barrel_diameter
        have_inner = config%cylinder%has_inner_diameter
        if (have_inner) inner_radius = 0.5_c_double * config%cylinder%inner_diameter
        have_axial_min = config%cylinder%has_axial_start
        if (have_axial_min) axial_min = config%cylinder%axial_start
        have_axial_max = have_axial_min .and. config%cylinder%has_barrel_length
        if (have_axial_max) axial_max = config%cylinder%axial_start + config%cylinder%barrel_length

        if (.not.(have_outer .or. have_inner .or. have_axial_min .or. have_axial_max)) return

        face_count = boundary_faces%count
        if (face_count <= 0) return
        allocate(face_is_outer(face_count))
        allocate(face_is_inner(face_count))
        allocate(face_is_axial_min(face_count))
        allocate(face_is_axial_max(face_count))
        face_is_outer = .false.
        face_is_inner = .false.
        face_is_axial_min = .false.
        face_is_axial_max = .false.

        do face_idx = 1, face_count
            call gather_face_vertices(boundary_faces%elements(face_idx), boundary_faces%face_ids(face_idx), kvert, node_ids)
            if (have_outer) then
                face_is_outer(face_idx) = check_all_vertices_on_radius(coords, node_ids, outer_radius, &
                    classification%tolerance)
            else
                face_is_outer(face_idx) = .false.
            end if
            if (have_inner) then
                face_is_inner(face_idx) = check_all_vertices_on_radius(coords, node_ids, inner_radius, &
                    classification%tolerance)
            else
                face_is_inner(face_idx) = .false.
            end if
            if (have_axial_min) then
                face_is_axial_min(face_idx) = check_all_vertices_on_plane(coords, node_ids, 3, axial_min, &
                    classification%tolerance)
            else
                face_is_axial_min(face_idx) = .false.
            end if
            if (have_axial_max) then
                face_is_axial_max(face_idx) = check_all_vertices_on_plane(coords, node_ids, 3, axial_max, &
                    classification%tolerance)
            else
                face_is_axial_max(face_idx) = .false.
            end if
        end do

        if (have_outer) then
            call subset_face_list(boundary_faces, face_is_outer, classification%cyl_outer_faces)
            call collect_nodes_from_face_list(classification%cyl_outer_faces, coords, kvert, classification%cyl_outer_nodes)
        end if
        if (have_inner) then
            call subset_face_list(boundary_faces, face_is_inner, classification%cyl_inner_faces)
            call collect_nodes_from_face_list(classification%cyl_inner_faces, coords, kvert, classification%cyl_inner_nodes)
        end if
        if (have_axial_min) then
            call subset_face_list(boundary_faces, face_is_axial_min, classification%axial_min_faces)
            call collect_nodes_from_face_list(classification%axial_min_faces, coords, kvert, classification%axial_min_nodes)
        end if
        if (have_axial_max) then
            call subset_face_list(boundary_faces, face_is_axial_max, classification%axial_max_faces)
            call collect_nodes_from_face_list(classification%axial_max_faces, coords, kvert, classification%axial_max_nodes)
        end if

        call subset_face_list(boundary_faces, .not.(face_is_outer .or. face_is_inner .or. face_is_axial_min .or. &
            face_is_axial_max), classification%all_boundary_faces)
        call collect_nodes_from_face_list(classification%all_boundary_faces, coords, kvert, classification%inner_wall_nodes)

        deallocate(face_is_outer, face_is_inner, face_is_axial_min, face_is_axial_max)
    end subroutine classify_hollowcylinder_boundaries

    subroutine classify_box_boundaries(config, coords, kvert, knpr, boundary_faces, classification)
        type(MeshConfig), intent(in) :: config
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: kvert(:, :)
        integer(c_int), intent(in) :: knpr(:)
        type(FaceList), intent(in) :: boundary_faces
        type(BoxBoundaryClassification), intent(out) :: classification

        real(c_double) :: start(3), lengths(3)
        real(c_double) :: plane_values(6)
        logical :: have_start, have_length
        integer :: face_count
        integer :: face_idx
        integer :: node_ids(4)
        logical, allocatable :: face_flags(:, :)

        classification%min_edge_length = compute_min_edge_length(coords, kvert)
        classification%tolerance = 0.5_c_double * classification%min_edge_length
        classification%total_boundary_nodes = count(knpr == 1_c_int)
        call allocate_zero_length_list(classification%x_min_nodes)
        call allocate_zero_length_list(classification%x_max_nodes)
        call allocate_zero_length_list(classification%y_min_nodes)
        call allocate_zero_length_list(classification%y_max_nodes)
        call allocate_zero_length_list(classification%z_min_nodes)
        call allocate_zero_length_list(classification%z_max_nodes)
        call allocate_zero_length_list(classification%inner_wall_nodes)
        call init_face_list(classification%x_min_faces)
        call init_face_list(classification%x_max_faces)
        call init_face_list(classification%y_min_faces)
        call init_face_list(classification%y_max_faces)
        call init_face_list(classification%z_min_faces)
        call init_face_list(classification%z_max_faces)
        call init_face_list(classification%all_boundary_faces)
        if (classification%total_boundary_nodes == 0) return

        have_start = config%box%has_geometry_start
        have_length = config%box%has_geometry_length
        if (.not. (have_start .and. have_length)) return

        face_count = boundary_faces%count
        if (face_count <= 0) return

        start = config%box%geometry_start
        lengths = config%box%geometry_length
        plane_values = [ &
            start(1), start(1) + lengths(1), &
            start(2), start(2) + lengths(2), &
            start(3), start(3) + lengths(3)]

        allocate(face_flags(6, face_count))
        face_flags = .false.

        do face_idx = 1, face_count
            call gather_face_vertices(boundary_faces%elements(face_idx), boundary_faces%face_ids(face_idx), kvert, node_ids)
            if (check_all_vertices_on_plane(coords, node_ids, 1, plane_values(1), classification%tolerance)) then
                face_flags(1, face_idx) = .true.
            end if
            if (check_all_vertices_on_plane(coords, node_ids, 1, plane_values(2), classification%tolerance)) then
                face_flags(2, face_idx) = .true.
            end if
            if (check_all_vertices_on_plane(coords, node_ids, 2, plane_values(3), classification%tolerance)) then
                face_flags(3, face_idx) = .true.
            end if
            if (check_all_vertices_on_plane(coords, node_ids, 2, plane_values(4), classification%tolerance)) then
                face_flags(4, face_idx) = .true.
            end if
            if (check_all_vertices_on_plane(coords, node_ids, 3, plane_values(5), classification%tolerance)) then
                face_flags(5, face_idx) = .true.
            end if
            if (check_all_vertices_on_plane(coords, node_ids, 3, plane_values(6), classification%tolerance)) then
                face_flags(6, face_idx) = .true.
            end if
        end do

        call subset_face_list(boundary_faces, face_flags(1, :), classification%x_min_faces)
        call subset_face_list(boundary_faces, face_flags(2, :), classification%x_max_faces)
        call subset_face_list(boundary_faces, face_flags(3, :), classification%y_min_faces)
        call subset_face_list(boundary_faces, face_flags(4, :), classification%y_max_faces)
        call subset_face_list(boundary_faces, face_flags(5, :), classification%z_min_faces)
        call subset_face_list(boundary_faces, face_flags(6, :), classification%z_max_faces)
        call collect_nodes_from_face_list(classification%x_min_faces, coords, kvert, classification%x_min_nodes)
        call collect_nodes_from_face_list(classification%x_max_faces, coords, kvert, classification%x_max_nodes)
        call collect_nodes_from_face_list(classification%y_min_faces, coords, kvert, classification%y_min_nodes)
        call collect_nodes_from_face_list(classification%y_max_faces, coords, kvert, classification%y_max_nodes)
        call collect_nodes_from_face_list(classification%z_min_faces, coords, kvert, classification%z_min_nodes)
        call collect_nodes_from_face_list(classification%z_max_faces, coords, kvert, classification%z_max_nodes)

        call subset_face_list(boundary_faces, .not.(face_flags(1, :) .or. face_flags(2, :) .or. face_flags(3, :) .or. &
            face_flags(4, :) .or. face_flags(5, :) .or. face_flags(6, :)), classification%all_boundary_faces)
        call collect_nodes_from_face_list(classification%all_boundary_faces, coords, kvert, classification%inner_wall_nodes)

        deallocate(face_flags)
    end subroutine classify_box_boundaries

    real(c_double) function compute_min_edge_length(coords, kvert) result(min_len)
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: kvert(:, :)
        integer :: nel, nve, elem_idx, edge_idx
        integer :: vid_a, vid_b
        real(c_double) :: diff(3), edge_len

        nel = size(kvert, 2)
        nve = size(kvert, 1)
        if (nel == 0 .or. nve < 2) then
            min_len = 0.0_c_double
            return
        end if

        min_len = huge(1.0_c_double)
        do elem_idx = 1, nel
            do edge_idx = 1, size(hexahedron_edges, 2)
                vid_a = int(kvert(hexahedron_edges(1, edge_idx), elem_idx))
                vid_b = int(kvert(hexahedron_edges(2, edge_idx), elem_idx))
                diff = coords(:, vid_a) - coords(:, vid_b)
                edge_len = sqrt(sum(diff * diff))
                if (edge_len > 0.0_c_double .and. edge_len < min_len) min_len = edge_len
            end do
        end do

        if (min_len == huge(1.0_c_double)) min_len = 0.0_c_double
    end function compute_min_edge_length

    subroutine allocate_zero_length_list(list)
        integer, allocatable, intent(inout) :: list(:)
        if (allocated(list)) deallocate(list)
        allocate(list(0))
    end subroutine allocate_zero_length_list

    subroutine init_face_list(list)
        type(FaceList), intent(inout) :: list
        if (allocated(list%elements)) deallocate(list%elements)
        if (allocated(list%face_ids)) deallocate(list%face_ids)
        list%count = 0
    end subroutine init_face_list

    subroutine allocate_face_list(list, count)
        type(FaceList), intent(inout) :: list
        integer, intent(in) :: count

        call init_face_list(list)
        if (count <= 0) then
            list%count = 0
            if (.not. allocated(list%elements)) allocate(list%elements(0))
            if (.not. allocated(list%face_ids)) allocate(list%face_ids(0))
            return
        end if
        list%count = count
        allocate(list%elements(count))
        allocate(list%face_ids(count))
    end subroutine allocate_face_list

    subroutine copy_face_list(src, dest)
        type(FaceList), intent(in) :: src
        type(FaceList), intent(out) :: dest

        call allocate_face_list(dest, src%count)
        if (src%count > 0) then
            dest%elements = src%elements
            dest%face_ids = src%face_ids
        end if
    end subroutine copy_face_list

    subroutine subset_face_list(src, flags, dest)
        type(FaceList), intent(in) :: src
        logical, intent(in) :: flags(:)
        type(FaceList), intent(out) :: dest
        integer :: cnt, idx, pos

        if (size(flags) /= src%count) then
            call allocate_face_list(dest, 0)
            return
        end if
        cnt = count(flags)
        call allocate_face_list(dest, cnt)
        if (cnt == 0) return
        pos = 0
        do idx = 1, src%count
            if (.not. flags(idx)) cycle
            pos = pos + 1
            dest%elements(pos) = src%elements(idx)
            dest%face_ids(pos) = src%face_ids(idx)
        end do
    end subroutine subset_face_list

    subroutine collect_nodes_from_mask(mask, out_list)
        logical, intent(in) :: mask(:)
        integer, allocatable, intent(out) :: out_list(:)
        integer :: cnt, idx, pos

        cnt = count(mask)
        if (allocated(out_list)) deallocate(out_list)
        allocate(out_list(cnt))
        if (cnt == 0) return
        pos = 0
        do idx = 1, size(mask)
            if (.not. mask(idx)) cycle
            pos = pos + 1
            out_list(pos) = idx
        end do
    end subroutine collect_nodes_from_mask

    subroutine gather_face_vertices(elem_idx, face_idx, kvert, node_ids)
        integer, intent(in) :: elem_idx, face_idx
        integer(c_int), intent(in) :: kvert(:, :)
        integer, intent(out) :: node_ids(4)
        integer :: i

        do i = 1, 4
            node_ids(i) = int(kvert(hexahedron_faces(i, face_idx), elem_idx))
        end do
    end subroutine gather_face_vertices

    subroutine compute_face_center(coords, node_ids, center)
        real(c_double), intent(in) :: coords(:, :)
        integer, intent(in) :: node_ids(4)
        real(c_double), intent(out) :: center(3)
        integer :: i

        center = 0.0_c_double
        do i = 1, 4
            center = center + coords(:, node_ids(i))
        end do
        center = center / 4.0_c_double
    end subroutine compute_face_center

    logical function check_all_vertices_on_plane(coords, node_ids, axis, plane_value, tolerance) result(all_on_plane)
        real(c_double), intent(in) :: coords(:, :)
        integer, intent(in) :: node_ids(4)
        integer, intent(in) :: axis
        real(c_double), intent(in) :: plane_value, tolerance
        integer :: i

        all_on_plane = .true.
        do i = 1, 4
            if (abs(coords(axis, node_ids(i)) - plane_value) > tolerance) then
                all_on_plane = .false.
                return
            end if
        end do
    end function check_all_vertices_on_plane

    logical function check_all_vertices_on_radius(coords, node_ids, target_radius, tolerance) result(all_on_radius)
        real(c_double), intent(in) :: coords(:, :)
        integer, intent(in) :: node_ids(4)
        real(c_double), intent(in) :: target_radius, tolerance
        integer :: i
        real(c_double) :: radial

        all_on_radius = .true.
        do i = 1, 4
            radial = sqrt(coords(1, node_ids(i))**2 + coords(2, node_ids(i))**2)
            if (abs(radial - target_radius) > tolerance) then
                all_on_radius = .false.
                return
            end if
        end do
    end function check_all_vertices_on_radius

    subroutine sort_face_vertices(values)
        integer, intent(inout) :: values(4)
        integer :: i, j, temp

        do i = 2, 4
            temp = values(i)
            j = i - 1
            do while (j >= 1 .and. values(j) > temp)
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = temp
        end do
    end subroutine sort_face_vertices

    subroutine sort_faces_by_keys(keys, order)
        integer, intent(in) :: keys(:, :)
        integer, intent(inout) :: order(:)

        call quicksort_faces(keys, order, 1, size(order))
    end subroutine sort_faces_by_keys

    recursive subroutine quicksort_faces(keys, order, left, right)
        integer, intent(in) :: keys(:, :)
        integer, intent(inout) :: order(:)
        integer, intent(in) :: left, right
        integer :: i, j, pivot_idx, temp

        if (left >= right) return
        i = left
        j = right
        pivot_idx = order((left + right) / 2)
        do
            do while (face_compare(keys, order(i), pivot_idx) < 0)
                i = i + 1
            end do
            do while (face_compare(keys, order(j), pivot_idx) > 0)
                j = j - 1
            end do
            if (i <= j) then
                temp = order(i)
                order(i) = order(j)
                order(j) = temp
                i = i + 1
                j = j - 1
            end if
            if (i > j) exit
        end do
        if (left < j) call quicksort_faces(keys, order, left, j)
        if (i < right) call quicksort_faces(keys, order, i, right)
    end subroutine quicksort_faces

    integer function face_compare(keys, idx_a, idx_b) result(cmp)
        integer, intent(in) :: keys(:, :)
        integer, intent(in) :: idx_a, idx_b
        integer :: k

        do k = 1, size(keys, 1)
            if (keys(k, idx_a) < keys(k, idx_b)) then
                cmp = -1
                return
            else if (keys(k, idx_a) > keys(k, idx_b)) then
                cmp = 1
                return
            end if
        end do
        cmp = 0
    end function face_compare

    logical function faces_equal(keys, idx_a, idx_b) result(is_equal)
        integer, intent(in) :: keys(:, :)
        integer, intent(in) :: idx_a, idx_b

        is_equal = all(keys(:, idx_a) == keys(:, idx_b))
    end function faces_equal

    subroutine collect_nodes_from_face_list(faces, coords, kvert, out_nodes)
        type(FaceList), intent(in) :: faces
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: kvert(:, :)
        integer, allocatable, intent(out) :: out_nodes(:)
        logical, allocatable :: mask(:)
        integer :: idx, node_ids(4)

        if (faces%count <= 0) then
            call allocate_zero_length_list(out_nodes)
            return
        end if

        allocate(mask(size(coords, 2)))
        mask = .false.
        do idx = 1, faces%count
            call gather_face_vertices(faces%elements(idx), faces%face_ids(idx), kvert, node_ids)
            mask(node_ids) = .true.
        end do
        call collect_nodes_from_mask(mask, out_nodes)
        deallocate(mask)
    end subroutine collect_nodes_from_face_list

end module bc_treatment
