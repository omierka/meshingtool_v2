program fortran_cgal_demo
    use iso_c_binding, only: c_ptr, c_null_char, c_associated, c_size_t, c_double, c_int, c_char
    implicit none

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

        function cgal_copy_vertices(handle, buffer, buffer_len) bind(C, name="cgal_copy_vertices") result(status)
            import :: c_ptr, c_size_t, c_double, c_int
            type(c_ptr), value :: handle
            real(c_double) :: buffer(*)
            integer(c_size_t), value :: buffer_len
            integer(c_int) :: status
        end function cgal_copy_vertices

        function cgal_copy_triangles(handle, buffer, buffer_len) bind(C, name="cgal_copy_triangles") result(status)
            import :: c_ptr, c_size_t, c_int
            type(c_ptr), value :: handle
            integer(c_int) :: buffer(*)
            integer(c_size_t), value :: buffer_len
            integer(c_int) :: status
        end function cgal_copy_triangles

        function cgal_minimal_bounding_sphere(points, point_count, center, radius) &
            bind(C, name="cgal_minimal_bounding_sphere") result(status)
            import :: c_size_t, c_double, c_int
            real(c_double) :: points(*)
            integer(c_size_t), value :: point_count
            real(c_double) :: center(*)
            real(c_double) :: radius
            integer(c_int) :: status
        end function cgal_minimal_bounding_sphere

        function cgal_sphere_intersects_mesh(handle, center, radius, intersects) &
            bind(C, name="cgal_sphere_intersects_mesh") result(status)
            import :: c_ptr, c_double, c_int
            type(c_ptr), value :: handle
            real(c_double) :: center(*)
            real(c_double), value :: radius
            integer(c_int) :: intersects
            integer(c_int) :: status
        end function cgal_sphere_intersects_mesh

        function cgal_point_inside_mesh(handle, point, inside) bind(C, name="cgal_point_inside_mesh") result(status)
            import :: c_ptr, c_double, c_int
            type(c_ptr), value :: handle
            real(c_double) :: point(*)
            integer(c_int) :: inside
            integer(c_int) :: status
        end function cgal_point_inside_mesh

        function cgal_closest_point_index_to_mesh(handle, points, point_count, closest_index) &
            bind(C, name="cgal_closest_point_index_to_mesh") result(status)
            import :: c_ptr, c_size_t, c_double, c_int
            type(c_ptr), value :: handle
            real(c_double) :: points(*)
            integer(c_size_t), value :: point_count
            integer(c_int) :: closest_index
            integer(c_int) :: status
        end function cgal_closest_point_index_to_mesh

        function cgal_triangle_intersects_mesh(handle, triangle_points, intersects) &
            bind(C, name="cgal_triangle_intersects_mesh") result(status)
            import :: c_ptr, c_double, c_int
            type(c_ptr), value :: handle
            real(c_double) :: triangle_points(*)
            integer(c_int) :: intersects
            integer(c_int) :: status
        end function cgal_triangle_intersects_mesh

        function cgal_segment_intersects_mesh(handle, segment_points, intersects) &
            bind(C, name="cgal_segment_intersects_mesh") result(status)
            import :: c_ptr, c_double, c_int
            type(c_ptr), value :: handle
            real(c_double) :: segment_points(*)
            integer(c_int) :: intersects
            integer(c_int) :: status
        end function cgal_segment_intersects_mesh
    end interface

    character(len=256) :: surface_file, hex_file, output_folder
    real(c_double) :: hex_scale
    character(len=256) :: header_line1, header_line2
    character(kind=c_char), allocatable :: c_path(:)
    type(c_ptr) :: mesh_handle
    integer(c_size_t) :: vertex_count, triangle_count
    real(c_double), allocatable :: vertices(:)
    integer(c_int), allocatable :: triangles(:)
    ! Hex mesh data
    integer :: hex_nel, hex_nvt, hex_nve
    integer :: hex_nbct, hex_nee, hex_nae
    real(c_double), allocatable :: hex_coords(:, :)
    integer(c_int), allocatable :: hex_kvert(:, :)
    integer(c_int), allocatable :: hex_knpr(:)
    integer(c_int) :: ierr
    integer(c_size_t) :: buffer_len
    real(c_double), allocatable :: elem_vertices(:)
    real(c_double), allocatable :: local_coords(:, :)
    real(c_double) :: sphere_center(3), sphere_radius
    real(c_double) :: vertex_point(3)
    integer(c_int) :: intersects_flag
    integer(c_int) :: inside_flag
    integer :: elem_idx, local_idx
    logical, allocatable :: intersect_mask(:)
    integer :: intersect_count
    real(c_double), allocatable :: filtered_coords(:, :)
    integer(c_int), allocatable :: filtered_kvert(:, :)
    integer(c_int), allocatable :: filtered_knpr(:)
    integer :: filtered_nel, filtered_nvt
    logical :: element_selected, vertex_inside, all_vertices_inside
    integer :: inside_vertex_count
    integer :: count_intersected, count_inside, count_edge_intersected, count_diagonal_intersected, count_discarded
    integer :: progress_pct, last_progress, progress_char
    character(len=100) :: progress_bar_template, progress_bar_line
    character(len=512) :: filtered_hex_file = "Filtered.tri"
    character(len=512) :: filtered_vtu_file = "Filtered.vtu"

    surface_file = "sphere.off"
    hex_file = "single.tri"
    hex_scale = 1.0_c_double
    output_folder= "."
    
    call parse_arguments(surface_file, hex_file, hex_scale, output_folder)

    c_path = to_c_string(trim(surface_file))
    mesh_handle = cgal_load_off(c_path)
    if (.not. c_associated(mesh_handle)) then
        error stop "CGAL failed to load OFF file."
    end if

    vertex_count = cgal_get_vertex_count(mesh_handle)
    triangle_count = cgal_get_triangle_count(mesh_handle)

    allocate(vertices(3 * vertex_count))
    allocate(triangles(3 * triangle_count))

    buffer_len = int(size(vertices), kind=c_size_t)
    ierr = cgal_copy_vertices(mesh_handle, vertices, buffer_len)
    if (ierr /= 0) then
        call cgal_free_mesh(mesh_handle)
        error stop "Copying vertices from CGAL failed."
    end if

    buffer_len = int(size(triangles), kind=c_size_t)
    ierr = cgal_copy_triangles(mesh_handle, triangles, buffer_len)
    if (ierr /= 0) then
        call cgal_free_mesh(mesh_handle)
        error stop "Copying triangle connectivity from CGAL failed."
    end if

    call read_single_tri(trim(hex_file), hex_nel, hex_nvt, hex_nve, hex_nbct, hex_nee, hex_nae, header_line1, header_line2, &
                         hex_coords, hex_kvert, hex_knpr)
    hex_coords = hex_scale * hex_coords

    allocate(elem_vertices(3 * hex_nve))
    allocate(local_coords(3, hex_nve))
    allocate(intersect_mask(hex_nel))
    intersect_mask = .false.
    count_intersected = 0
    count_inside = 0
    count_edge_intersected = 0
    count_diagonal_intersected = 0
    last_progress = 0
    progress_bar_template = ''
    progress_bar_line = ''
    do progress_char = 1, len(progress_bar_template)
        progress_bar_template(progress_char:progress_char) = '%'
        progress_bar_line(progress_char:progress_char) = ' '
    end do
    write (*, '(A)') progress_bar_template
    write (*, '(A)', advance='no') progress_bar_line
    do elem_idx = 1, hex_nel
        call gather_element_vertices(elem_idx, elem_vertices)
        do local_idx = 1, hex_nve
            local_coords(1, local_idx) = elem_vertices(3 * (local_idx - 1) + 1)
            local_coords(2, local_idx) = elem_vertices(3 * (local_idx - 1) + 2)
            local_coords(3, local_idx) = elem_vertices(3 * (local_idx - 1) + 3)
        end do
        ierr = cgal_minimal_bounding_sphere(elem_vertices, int(hex_nve, kind=c_size_t), sphere_center, sphere_radius)
        if (ierr /= 0) then
            call cgal_free_mesh(mesh_handle)
            write (*, *) "Failed sphere computation for element ", elem_idx
            error stop "CGAL minimal bounding sphere computation failed."
        end if
        ierr = cgal_sphere_intersects_mesh(mesh_handle, sphere_center, sphere_radius, intersects_flag)
        if (ierr /= 0) then
            call cgal_free_mesh(mesh_handle)
            write (*, *) "Failed intersection check for element ", elem_idx
            error stop "CGAL sphere/mesh intersection check failed."
        end if
        inside_vertex_count = 0
        do local_idx = 1, hex_nve
            vertex_point = local_coords(:, local_idx)
            ierr = cgal_point_inside_mesh(mesh_handle, vertex_point, inside_flag)
            if (ierr /= 0) then
                call cgal_free_mesh(mesh_handle)
                write (*, *) "Failed vertex inside check for element ", elem_idx
                error stop "CGAL point-inside-mesh check failed."
            end if
            if (inside_flag == 1) then
                inside_vertex_count = inside_vertex_count + 1
            end if
        end do
        vertex_inside = (inside_vertex_count > 0)
        all_vertices_inside = (inside_vertex_count == hex_nve)

        element_selected = .false.
        if (all_vertices_inside) then
            element_selected = .true.
            count_inside = count_inside + 1
        else if (vertex_inside) then
            element_selected = .true.
            count_intersected = count_intersected + 1
        else
            if (intersects_flag == 1) then
                call filter_intersection_via_edges(mesh_handle, elem_idx, local_coords, element_selected)
                if (element_selected) then
                    count_edge_intersected = count_edge_intersected + 1
                else
                    call filter_intersection_via_diagonals(mesh_handle, elem_idx, local_coords, element_selected)
                    if (element_selected) count_diagonal_intersected = count_diagonal_intersected + 1
                end if
            else
                ierr = cgal_point_inside_mesh(mesh_handle, sphere_center, inside_flag)
                if (ierr /= 0) then
                    call cgal_free_mesh(mesh_handle)
                    write (*, *) "Failed sphere-center inside check for element ", elem_idx
                    error stop "CGAL point-inside-mesh check failed."
                end if
                if (inside_flag == 1) then
                    element_selected = .true.
                    count_inside = count_inside + 1
                end if
            end if
        end if
        if (element_selected) intersect_mask(elem_idx) = .true.

        progress_pct = int(100.0_c_double * real(elem_idx, kind=c_double) / real(hex_nel, kind=c_double))
        if (progress_pct > 100) progress_pct = 100
        if (progress_pct > last_progress) then
            do progress_char = max(1, last_progress + 1), progress_pct
                if (progress_char > len(progress_bar_line)) exit
                progress_bar_line(progress_char:progress_char) = '%'
            end do
            write (*, '(A)', advance='no') char(13)//progress_bar_line
            last_progress = progress_pct
        end if
    end do
    intersect_count = count(intersect_mask)
    count_discarded = hex_nel - (count_inside + count_intersected + count_edge_intersected + count_diagonal_intersected)

    if (intersect_count > 0) then
        write(filtered_hex_file(1:),'(A)') adjustl(trim(output_folder))//"/Filtered.tri"
        write(filtered_vtu_file(1:),'(A)') adjustl(trim(output_folder))//"/Filtered.vtu"
        
        call reduce_hex_mesh(intersect_mask, hex_coords, hex_kvert, hex_knpr, hex_nve, filtered_coords, filtered_kvert, &
                             filtered_knpr, filtered_nel, filtered_nvt)
        call write_single_tri(filtered_hex_file, header_line1, header_line2, filtered_nel, filtered_nvt, hex_nbct, &
                              hex_nve, hex_nee, hex_nae, filtered_coords, filtered_kvert, filtered_knpr)
        call write_vtu(filtered_vtu_file, filtered_coords, filtered_kvert, filtered_knpr)
    end if

    call cgal_free_mesh(mesh_handle)
    write (*,*)

    write (*, '(A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0)') '[DISCARDED]=', count_discarded, '[INTERSECTED]=', &
        count_intersected, '[INSIDE]=', count_inside, '[EDGEINTERSECTED]=', count_edge_intersected, &
        '[DIAGONALINTERSECTION]=', count_diagonal_intersected

contains
    subroutine parse_arguments(surface_file, hex_file, hex_scale, output_folder )
        character(len=*), intent(inout) :: surface_file, hex_file, output_folder
        real(c_double), intent(inout) :: hex_scale
        integer :: argc, i
        character(len=256) :: arg, value

        argc = command_argument_count()
        i = 1
        do while (i <= argc)
            call get_command_argument(i, arg)
            select case (trim(arg))
            case ("-t", "--tri")
                if (i == argc) call argument_error("Missing value for -t/--tri.")
                call get_command_argument(i + 1, value)
                surface_file = trim(adjustl(value))
                i = i + 2
            case ("-h", "--hex")
                if (i == argc) call argument_error("Missing value for -h/--hex.")
                call get_command_argument(i + 1, value)
                hex_file = trim(adjustl(value))
                i = i + 2
            case ("-s", "--hex-scale")
                if (i == argc) call argument_error("Missing value for -s/--hex-scale.")
                call get_command_argument(i + 1, value)
                read(value, *) hex_scale
                i = i + 2
           case ("-o", "--output-folder")
                if (i == argc) call argument_error("Missing value for -o/--output-folder.")
                call get_command_argument(i + 1, value)
                read(value, *)  output_folder
                i = i + 2
            case default
                call argument_error("Unknown argument: "//trim(arg))
            end select
        end do
    end subroutine parse_arguments

    subroutine argument_error(message)
        character(len=*), intent(in) :: message
        write (*, *) trim(message)
        call print_usage()
        error stop
    end subroutine argument_error

    subroutine print_usage()
        write (*, *) "Usage: fortran_cgal [-h hex_mesh.tri] [-t surface.off]"
    end subroutine print_usage

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

    subroutine read_single_tri(filename, nel, nvt, nve, nbct, nee, nae, header1, header2, coords, kvert, knpr)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: nel, nvt, nve, nbct, nee, nae
        character(len=*), intent(out) :: header1, header2
        real(c_double), allocatable, intent(out) :: coords(:, :)
        integer(c_int), allocatable, intent(out) :: kvert(:, :)
        integer(c_int), allocatable, intent(out) :: knpr(:)

        integer :: unit, i
        character(len=256) :: label

        open(newunit=unit, file=filename, status="old", action="read")
        read(unit, '(A)') header1
        read(unit, '(A)') header2
        read(unit, *) nel, nvt, nbct, nve, nee, nae

        read(unit, '(A)') label
        if (index(adjustl(label), "DCORVG") /= 1) error stop "Expected DCORVG section"
        allocate(coords(3, nvt))
        do i = 1, nvt
            read(unit, *) coords(:, i)
        end do

        read(unit, '(A)') label
        if (index(adjustl(label), "KVERT") /= 1) error stop "Expected KVERT section"
        allocate(kvert(nve, nel))
        do i = 1, nel
            read(unit, *) kvert(:, i)
        end do

        read(unit, '(A)') label
        if (index(adjustl(label), "KNPR") /= 1) error stop "Expected KNPR section"
        allocate(knpr(nvt))
        do i = 1, nvt
            read(unit, *) knpr(i)
        end do

        close(unit)
    end subroutine read_single_tri

    subroutine reduce_hex_mesh(intersect_mask, coords, kvert, knpr, nve, out_coords, out_kvert, out_knpr, out_nel, out_nvt)
        logical, intent(in) :: intersect_mask(:)
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: kvert(:, :)
        integer(c_int), intent(in) :: knpr(:)
        integer, intent(in) :: nve
        real(c_double), allocatable, intent(out) :: out_coords(:, :)
        integer(c_int), allocatable, intent(out) :: out_kvert(:, :)
        integer(c_int), allocatable, intent(out) :: out_knpr(:)
        integer, intent(out) :: out_nel, out_nvt

        integer :: hex_nel, hex_nvt
        logical, allocatable :: vertex_used(:)
        integer, allocatable :: vertex_remap(:)
        integer, allocatable :: temp_kvert(:, :)
        integer :: elem_idx, local_idx, new_elem_idx, vid, new_vid

        hex_nel = size(kvert, 2)
        hex_nvt = size(coords, 2)
        out_nel = count(intersect_mask)
        if (out_nel == 0) then
            out_nvt = 0
            allocate(out_coords(3, 0))
            allocate(out_kvert(nve, 0))
            allocate(out_knpr(0))
            return
        end if

        allocate(vertex_used(hex_nvt))
        vertex_used = .false.
        allocate(temp_kvert(nve, out_nel))

        new_elem_idx = 0
        do elem_idx = 1, hex_nel
            if (intersect_mask(elem_idx)) then
                new_elem_idx = new_elem_idx + 1
                do local_idx = 1, nve
                    vid = int(kvert(local_idx, elem_idx))
                    vertex_used(vid) = .true.
                    temp_kvert(local_idx, new_elem_idx) = vid
                end do
            end if
        end do

        allocate(vertex_remap(hex_nvt))
        vertex_remap = 0
        out_nvt = 0
        do vid = 1, hex_nvt
            if (vertex_used(vid)) then
                out_nvt = out_nvt + 1
                vertex_remap(vid) = out_nvt
            end if
        end do

        allocate(out_coords(3, out_nvt))
        allocate(out_knpr(out_nvt))
        do vid = 1, hex_nvt
            new_vid = vertex_remap(vid)
            if (new_vid > 0) then
                out_coords(:, new_vid) = coords(:, vid)
                out_knpr(new_vid) = knpr(vid)
            end if
        end do

        allocate(out_kvert(nve, out_nel))
        do elem_idx = 1, out_nel
            do local_idx = 1, nve
                vid = temp_kvert(local_idx, elem_idx)
                out_kvert(local_idx, elem_idx) = int(vertex_remap(vid), kind=c_int)
            end do
        end do

        deallocate(vertex_used, vertex_remap, temp_kvert)
    end subroutine reduce_hex_mesh

    subroutine write_single_tri(filename, header1, header2, nel, nvt, nbct, nve, nee, nae, coords, kvert, knpr)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: header1, header2
        integer, intent(in) :: nel, nvt, nbct, nve, nee, nae
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: kvert(:, :)
        integer(c_int), intent(in) :: knpr(:)

        integer :: unit, i

        open(newunit=unit, file=filename, status="replace", action="write")
        write(unit, '(A)') trim(header1)
        write(unit, '(A)') trim(header2)
        write(unit, '(I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,3X,A)') nel, nvt, nbct, nve, nee, nae, "NEL,NVT,NBCT,NVE,NEE,NAE"

        write(unit, '(A)') "DCORVG"
        do i = 1, nvt
            write(unit, '(3(1X,ES24.16))') coords(1, i), coords(2, i), coords(3, i)
        end do

        write(unit, '(A)') "KVERT"
        do i = 1, nel
            write(unit, '(*(1X,I10))') kvert(:, i)
        end do

        write(unit, '(A)') "KNPR"
        do i = 1, nvt
            write(unit, '(I2)') knpr(i)
        end do

        close(unit)
    end subroutine write_single_tri

    subroutine write_vtu(filename, coords, kvert, knpr)
        character(len=*), intent(in) :: filename
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: kvert(:, :)
        integer(c_int), intent(in) :: knpr(:)

        integer :: unit, nvt, nel, nve
        integer :: i, j
        integer :: offset

        nvt = size(coords, 2)
        nel = size(kvert, 2)
        nve = size(kvert, 1)

        open(newunit=unit, file=filename, status="replace", action="write")
        write(unit, '(A)') '<?xml version="1.0"?>'
        write(unit, '(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
        write(unit, '(A)') '  <UnstructuredGrid>'
        write(unit, '(A,I0,A,I0,A)') '    <Piece NumberOfPoints="', nvt, '" NumberOfCells="', nel, '">'

        write(unit, '(A)') '      <PointData Scalars="KNPR">'
        write(unit, '(A)') '        <DataArray type="Int32" Name="KNPR" format="ascii">'
        write(unit, '(6(1X,I8))') (int(knpr(i), kind=4), i = 1, nvt)
        write(unit, '(A)') '        </DataArray>'
        write(unit, '(A)') '      </PointData>'

        write(unit, '(A)') '      <CellData/>'

        write(unit, '(A)') '      <Points>'
        write(unit, '(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
        do i = 1, nvt
            write(unit, '(3(1X,ES24.16))') coords(1, i), coords(2, i), coords(3, i)
        end do
        write(unit, '(A)') '        </DataArray>'
        write(unit, '(A)') '      </Points>'

        write(unit, '(A)') '      <Cells>'
        write(unit, '(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
        do i = 1, nel
            write(unit, '(*(1X,I10))') (int(kvert(j, i) - 1, kind=4), j = 1, nve)
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
    end subroutine write_vtu

    subroutine filter_intersection_via_edges(mesh_handle, elem_idx, coords, element_selected)
        type(c_ptr), intent(in) :: mesh_handle
        integer, intent(in) :: elem_idx
        real(c_double), intent(in) :: coords(:, :)
        logical, intent(out) :: element_selected

        integer, parameter :: hexahedron_edges(2, 12) = reshape([ &
            1, 2, 2, 3, 3, 4, 4, 1, 5, 6, 6, 7, 7, 8, 8, 5, 1, 5, 2, 6, 3, 7, 4, 8], [2, 12])
        integer :: edge_idx, v_start, v_end, nverts
        real(c_double) :: segment_points(6)
        integer(c_int) :: status_c, intersects_flag

        element_selected = .false.
        nverts = size(coords, 2)
        if (nverts < 8) return

        do edge_idx = 1, size(hexahedron_edges, 2)
            v_start = hexahedron_edges(1, edge_idx)
            v_end = hexahedron_edges(2, edge_idx)
            if (v_start > nverts .or. v_end > nverts) cycle

            segment_points(1:3) = coords(:, v_start)
            segment_points(4:6) = coords(:, v_end)
            status_c = cgal_segment_intersects_mesh(mesh_handle, segment_points, intersects_flag)
            if (status_c /= 0) then
                call cgal_free_mesh(mesh_handle)
                write (*, *) "Segment intersection check failed for element ", elem_idx
                error stop "CGAL segment intersection failed."
            end if
            if (intersects_flag == 1) then
                element_selected = .true.
                return
            end if
        end do
    end subroutine filter_intersection_via_edges

    subroutine filter_intersection_via_diagonals(mesh_handle, elem_idx, coords, element_selected)
        type(c_ptr), intent(in) :: mesh_handle
        integer, intent(in) :: elem_idx
        real(c_double), intent(in) :: coords(:, :)
        logical, intent(out) :: element_selected

        integer, parameter :: diagonal_pairs(2, 4) = reshape([1, 7, 2, 8, 3, 5, 4, 6], [2, 4])
        integer :: pair_idx, v_start, v_end, nverts
        real(c_double) :: segment_points(6)
        integer(c_int) :: status_c, intersects_flag

        element_selected = .false.
        nverts = size(coords, 2)
        if (nverts < 8) return

        do pair_idx = 1, size(diagonal_pairs, 2)
            v_start = diagonal_pairs(1, pair_idx)
            v_end = diagonal_pairs(2, pair_idx)
            if (v_start > nverts .or. v_end > nverts) cycle

            segment_points(1:3) = coords(:, v_start)
            segment_points(4:6) = coords(:, v_end)
            status_c = cgal_segment_intersects_mesh(mesh_handle, segment_points, intersects_flag)
            if (status_c /= 0) then
                call cgal_free_mesh(mesh_handle)
                write (*, *) "Diagonal intersection check failed for element ", elem_idx
                error stop "CGAL segment intersection failed."
            end if
            if (intersects_flag == 1) then
                element_selected = .true.
                return
            end if
        end do
    end subroutine filter_intersection_via_diagonals

    subroutine filter_intersection_via_faces(mesh_handle, elem_idx, ref_idx_c, coords, element_selected)
        type(c_ptr), intent(in) :: mesh_handle
        integer, intent(in) :: elem_idx
        integer(c_int), intent(in) :: ref_idx_c
        real(c_double), intent(in) :: coords(:, :)
        logical, intent(out) :: element_selected

        real(c_double) :: best_dists(3)
        integer :: neighbor_idx(3)
        real(c_double) :: dist2
        integer :: ref_idx, idx, insert_pos, shift_idx
        integer, dimension(2, 3) :: pair_map
        integer :: pair_first, pair_second, pair
        real(c_double) :: triangle_points(9)
        integer(c_int) :: status_c, intersects_flag

        element_selected = .false.
        ref_idx = int(ref_idx_c)
        if (ref_idx < 1 .or. ref_idx > hex_nve) return

        best_dists = huge(1.0_c_double)
        neighbor_idx = 0
        do idx = 1, hex_nve
            if (idx == ref_idx) cycle
            dist2 = sum((coords(:, idx) - coords(:, ref_idx))**2)
            insert_pos = 0
            do shift_idx = 1, 3
                if (dist2 < best_dists(shift_idx)) then
                    insert_pos = shift_idx
                    exit
                end if
            end do
            if (insert_pos > 0) then
                do shift_idx = 3, insert_pos + 1, -1
                    best_dists(shift_idx) = best_dists(shift_idx - 1)
                    neighbor_idx(shift_idx) = neighbor_idx(shift_idx - 1)
                end do
                best_dists(insert_pos) = dist2
                neighbor_idx(insert_pos) = idx
            end if
        end do

        if (any(neighbor_idx == 0)) return

        pair_map = reshape([1, 2, 1, 3, 2, 3], [2, 3])
        do pair = 1, 3
            pair_first = neighbor_idx(pair_map(1, pair))
            pair_second = neighbor_idx(pair_map(2, pair))
            triangle_points(1:3) = coords(:, ref_idx)
            triangle_points(4:6) = coords(:, pair_first)
            triangle_points(7:9) = coords(:, pair_second)
            status_c = cgal_triangle_intersects_mesh(mesh_handle, triangle_points, intersects_flag)
            if (status_c /= 0) then
                call cgal_free_mesh(mesh_handle)
                write (*, *) "Triangle intersection check failed for element ", elem_idx
                error stop "CGAL triangle intersection failed."
            end if
            if (intersects_flag == 1) then
                element_selected = .true.
                return
            end if
        end do
    end subroutine filter_intersection_via_faces

    subroutine gather_element_vertices(elem_idx, buffer)
        integer, intent(in) :: elem_idx
        real(c_double), intent(out) :: buffer(:)
        integer :: local_idx, base_idx, vertex_id

        if (size(buffer) /= 3 * hex_nve) error stop "gather_element_vertices buffer mismatch."
        if (elem_idx < 1 .or. elem_idx > hex_nel) error stop "Element index out of range."

        do local_idx = 1, hex_nve
            vertex_id = int(hex_kvert(local_idx, elem_idx))
            base_idx = 3 * (local_idx - 1)
            buffer(base_idx + 1) = hex_coords(1, vertex_id)
            buffer(base_idx + 2) = hex_coords(2, vertex_id)
            buffer(base_idx + 3) = hex_coords(3, vertex_id)
        end do
    end subroutine gather_element_vertices
end program fortran_cgal_demo
