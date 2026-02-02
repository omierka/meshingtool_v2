program fortran_cgal_demo
    use iso_c_binding, only: c_ptr, c_null_char, c_associated, c_size_t, c_double, c_int, c_char
    use iso_fortran_env, only: output_unit
    use mpi
    use bc_treatment, only: MeshConfig, initialize_mesh_config, load_mesh_config, log_mesh_config, &
        recompute_knpr_from_connectivity, HollowCylinderBoundaryClassification, classify_hollowcylinder_boundaries, &
        BoxBoundaryClassification, classify_box_boundaries, FaceList
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
    character(len=512) :: config_file
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
    integer :: local_count_intersected, local_count_inside, local_count_edge_intersected, local_count_diagonal_intersected
    integer :: worker_total_elems, worker_processed, worker_pct_sent
    real(c_double) :: start_time, elapsed_time
    character(len=512) :: filtered_tri_file = "Filtered.tri"
    character(len=512) :: filtered_hex_file = "Filtered.tri"
    character(len=512) :: filtered_vtu_file = "Filtered.vtu"
    type(MeshConfig) :: mesh_config
    type(HollowCylinderBoundaryClassification) :: hc_boundary
    type(BoxBoundaryClassification) :: box_boundary
    type(FaceList) :: boundary_faces
    logical :: hc_summary_ready, box_summary_ready
    character(len=512) :: mesh_output_folder
    character(len=512), allocatable :: meshdir_files(:)
    integer :: meshdir_file_count
    integer :: mpi_rank, mpi_size, mpi_err
    integer :: worker_count, worker_rank
    logical :: is_master, is_worker
    integer, parameter :: PROGRESS_TAG = 9001

    surface_file = "sphere.off"
    hex_file = "single.tri"
    hex_scale = 1.0_c_double
    output_folder= "."
    config_file = ""
    call initialize_mesh_config(mesh_config)
    hc_summary_ready = .false.
    box_summary_ready = .false.
    meshdir_file_count = 0
    allocate(meshdir_files(0))

    call MPI_Init(mpi_err)
    if (mpi_err /= MPI_SUCCESS) then
        error stop "MPI initialization failed."
    end if
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)
    call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, mpi_err)
    is_master = (mpi_rank == 0)
    worker_count = mpi_size - 1
    if (worker_count <= 0) then
        if (is_master) write (*,*) "This program requires MPI runs with at least 2 ranks."
        call MPI_Finalize(mpi_err)
        stop 1
    end if
    is_worker = (.not. is_master)
    worker_rank = mpi_rank - 1

    call parse_arguments(surface_file, hex_file, hex_scale, output_folder)
    call build_config_path(output_folder, "setup.e3d", config_file)
    call build_config_path(output_folder, "meshDir", mesh_output_folder)
    if (len_trim(config_file) > 0) then
        call load_mesh_config(trim(config_file), mesh_config)
        if (is_master) call log_mesh_config(mesh_config)
    end if

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
    local_count_intersected = 0
    local_count_inside = 0
    local_count_edge_intersected = 0
    local_count_diagonal_intersected = 0
    if (is_master) then
        write (*, '(A,I0,A)') 'Distributing filtering across ', mpi_size, ' MPI ranks (rank 0 orchestrates only).'
        write (*,'(A)', advance='no') 'Progress: '
        call flush(output_unit)
    end if
    start_time = MPI_Wtime()
    if (is_worker) then
        worker_total_elems = 0
        do elem_idx = worker_rank + 1, hex_nel, worker_count
            worker_total_elems = worker_total_elems + 1
        end do
        if (worker_total_elems <= 0) worker_total_elems = 1
        worker_processed = 0
        worker_pct_sent = 0
        do elem_idx = worker_rank + 1, hex_nel, worker_count
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
                local_count_inside = local_count_inside + 1
            else if (vertex_inside) then
                element_selected = .true.
                local_count_intersected = local_count_intersected + 1
            else
                if (intersects_flag == 1) then
                    call filter_intersection_via_edges(mesh_handle, elem_idx, local_coords, element_selected)
                    if (element_selected) then
                        local_count_edge_intersected = local_count_edge_intersected + 1
                    else
                        call filter_intersection_via_diagonals(mesh_handle, elem_idx, local_coords, element_selected)
                        if (element_selected) local_count_diagonal_intersected = local_count_diagonal_intersected + 1
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
                        local_count_inside = local_count_inside + 1
                    end if
                end if
            end if
            if (element_selected) intersect_mask(elem_idx) = .true.
            worker_processed = worker_processed + 1
            call report_worker_progress(worker_total_elems, worker_processed, worker_pct_sent)
        end do
        do while (worker_pct_sent < 100)
            worker_pct_sent = worker_pct_sent + 1
            call send_worker_progress(worker_pct_sent)
        end do
    end if
    if (is_master) call monitor_worker_progress(worker_count)
    if (mpi_size > 1) then
        if (is_master) then
            call MPI_Reduce(MPI_IN_PLACE, intersect_mask, hex_nel, MPI_LOGICAL, MPI_LOR, 0, MPI_COMM_WORLD, mpi_err)
        else
            call MPI_Reduce(intersect_mask, intersect_mask, hex_nel, MPI_LOGICAL, MPI_LOR, 0, MPI_COMM_WORLD, mpi_err)
        end if
        call MPI_Reduce(local_count_inside, count_inside, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_Reduce(local_count_intersected, count_intersected, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_Reduce(local_count_edge_intersected, count_edge_intersected, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_Reduce(local_count_diagonal_intersected, count_diagonal_intersected, 1, MPI_INTEGER, MPI_SUM, 0, &
            MPI_COMM_WORLD, mpi_err)
    else
        count_inside = local_count_inside
        count_intersected = local_count_intersected
        count_edge_intersected = local_count_edge_intersected
        count_diagonal_intersected = local_count_diagonal_intersected
    end if

    if (is_master) then
        intersect_count = count(intersect_mask)
        count_discarded = hex_nel - (count_inside + count_intersected + count_edge_intersected + count_diagonal_intersected)
    end if

    if (is_master .and. intersect_count > 0) then
        call ensure_directory_exists(mesh_output_folder)
        call build_config_path(trim(mesh_output_folder), "Filtered.tri", filtered_hex_file)
        write(filtered_vtu_file(1:),'(A)') adjustl(trim(output_folder))//"/Filtered.vtu"
        write(filtered_tri_file(1:),'(A)') adjustl(trim(output_folder))//"/Filtered.tri"
        call reduce_hex_mesh(intersect_mask, hex_coords, hex_kvert, hex_knpr, hex_nve, filtered_coords, filtered_kvert, &
                             filtered_knpr, filtered_nel, filtered_nvt)
        call recompute_knpr_from_connectivity(filtered_kvert, filtered_knpr, boundary_faces)
        if (mesh_config%loaded) then
            select case (trim(mesh_config%mesh_type))
            case ("HollowCylinder")
                call classify_hollowcylinder_boundaries(mesh_config, filtered_coords, filtered_kvert, filtered_knpr, &
                    boundary_faces, hc_boundary)
                hc_summary_ready = .true.
                call write_hc_parametrizations(mesh_output_folder, hc_boundary, meshdir_files, meshdir_file_count)
            case ("Box")
                call classify_box_boundaries(mesh_config, filtered_coords, filtered_kvert, filtered_knpr, boundary_faces, &
                    box_boundary)
                box_summary_ready = .true.
                call write_box_parametrizations(mesh_output_folder, box_boundary, meshdir_files, meshdir_file_count)
            end select
        end if
        call write_single_tri(filtered_tri_file, header_line1, header_line2, filtered_nel, filtered_nvt, hex_nbct, &
                              hex_nve, hex_nee, hex_nae, 1.0_c_double * filtered_coords, filtered_kvert, filtered_knpr)
        call write_single_tri(filtered_hex_file, header_line1, header_line2, filtered_nel, filtered_nvt, hex_nbct, &
                              hex_nve, hex_nee, hex_nae, 0.1_c_double * filtered_coords, filtered_kvert, filtered_knpr)
        call append_file_record(meshdir_files, meshdir_file_count, "Filtered.tri")
        call write_vtu(filtered_vtu_file, filtered_coords, filtered_kvert, filtered_knpr)
        call write_project_file(mesh_output_folder, meshdir_files, meshdir_file_count)
    end if

    if (is_master) elapsed_time = MPI_Wtime() - start_time

    call cgal_free_mesh(mesh_handle)

    if (is_master) then
        write (*, '(A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0)') '[DISCARDED]=', count_discarded, '[INTERSECTED]=', &
            count_intersected, '[INSIDE]=', count_inside, '[EDGEINTERSECTED]=', count_edge_intersected, &
            '[DIAGONALINTERSECTION]=', count_diagonal_intersected
        if (hc_summary_ready) call report_hollowcylinder_summary(hc_boundary)
        if (box_summary_ready) call report_box_summary(box_boundary)
        if (mpi_size > 1) then
            write (*,'(A,I0,A,F8.2,A)') 'MPI filtering completed across ', mpi_size, ' ranks in ', elapsed_time, ' s.'
        end if
    end if

    call MPI_Finalize(mpi_err)

contains
    subroutine report_worker_progress(total_elems, processed_elems, pct_sent)
        integer, intent(in) :: total_elems, processed_elems
        integer, intent(inout) :: pct_sent
        integer :: target_pct

        target_pct = int(100.0_c_double * real(processed_elems, kind=c_double) / &
            real(max(1, total_elems), kind=c_double))
        if (target_pct > 100) target_pct = 100
        do while (pct_sent < target_pct)
            pct_sent = pct_sent + 1
            call send_worker_progress(pct_sent)
        end do
    end subroutine report_worker_progress

    subroutine send_worker_progress(pct_value)
        integer, intent(in) :: pct_value
        integer :: ierr

        call MPI_Send(pct_value, 1, MPI_INTEGER, 0, PROGRESS_TAG, MPI_COMM_WORLD, ierr)
    end subroutine send_worker_progress

    subroutine monitor_worker_progress(worker_count)
        integer, intent(in) :: worker_count
        integer :: expected_messages, received_messages
        integer :: pct_value, pct_loop, status(MPI_STATUS_SIZE)
        integer :: global_pct, last_pct

        if (worker_count <= 0) then
            write (*,*)
            return
        end if

        expected_messages = worker_count * 100
        last_pct = 0
        global_pct = 0

        do received_messages = 1, expected_messages
            call MPI_Recv(pct_value, 1, MPI_INTEGER, MPI_ANY_SOURCE, PROGRESS_TAG, MPI_COMM_WORLD, status, mpi_err)
            global_pct = received_messages / worker_count
            if (global_pct > 100) global_pct = 100
            if (global_pct > last_pct) then
                do pct_loop = last_pct + 1, global_pct
                    write (*,'(A)', advance='no') '%'
                end do
                call flush(output_unit)
                last_pct = global_pct
            end if
        end do
        write (*,*)
    end subroutine monitor_worker_progress

    subroutine parse_arguments(surface_file, hex_file, hex_scale, output_folder)
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
                output_folder = trim(adjustl(value))
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

    subroutine build_config_path(folder, filename, out_path)
        character(len=*), intent(in) :: folder, filename
        character(len=*), intent(out) :: out_path
        character(len=len(out_path)) :: temp
        integer :: last

        temp = ""
        if (len_trim(folder) == 0) then
            temp = trim(filename)
        else
            temp = trim(adjustl(folder))
            last = len_trim(temp)
            if (last > 0) then
                if (temp(last:last) == '/' .or. temp(last:last) == '\') then
                    temp = trim(temp)//trim(filename)
                else
                    temp = trim(temp)//"/"//trim(filename)
                end if
            else
                temp = trim(filename)
            end if
        end if
        out_path = temp
    end subroutine build_config_path

    subroutine ensure_directory_exists(path)
        character(len=*), intent(in) :: path
        logical :: exists
        integer :: status
        character(len=1024) :: command

        inquire(file=trim(path), exist=exists)
        if (exists) return
        command = 'mkdir -p "'//trim(path)//'"'
        call execute_command_line(trim(command), exitstat=status)
        if (status /= 0) then
            write(*,'(A)') 'Warning: failed to create directory '//trim(path)
        end if
    end subroutine ensure_directory_exists

    subroutine append_file_record(records, record_count, filepath)
        character(len=512), allocatable, intent(inout) :: records(:)
        integer, intent(inout) :: record_count
        character(len=*), intent(in) :: filepath
        character(len=512), allocatable :: new_records(:)
        integer :: i

        allocate(new_records(record_count + 1))
        do i = 1, record_count
            new_records(i) = records(i)
        end do
        new_records(record_count + 1) = trim(filepath)
        record_count = record_count + 1
        call move_alloc(new_records, records)
    end subroutine append_file_record

    subroutine write_project_file(folder, records, record_count)
        character(len=*), intent(in) :: folder
        character(len=512), allocatable, intent(in) :: records(:)
        integer, intent(in) :: record_count
        character(len=512) :: filepath
        integer :: unit, i

        call build_config_path(folder, "file.prj", filepath)
        open(newunit=unit, file=trim(filepath), status="replace", action="write")
        do i = 1, record_count
            write(unit,'(A)') trim(records(i))
        end do
        close(unit)
    end subroutine write_project_file

    subroutine report_hollowcylinder_summary(classification)
        type(HollowCylinderBoundaryClassification), intent(in) :: classification
        integer :: count_outer, count_inner, count_zmin, count_zmax

        count_outer = size(classification%cyl_outer_nodes)
        count_inner = size(classification%cyl_inner_nodes)
        count_zmin = size(classification%axial_min_nodes)
        count_zmax = size(classification%axial_max_nodes)

        write(*,'("HC boundary tolerance=",ES12.5)') classification%tolerance
        write(*,'("HC boundary nodes: cyl_out=",I0," cyl_in=",I0," z-=",I0," z+=",I0," innerwall=",I0)') &
            count_outer, count_inner, count_zmin, count_zmax, size(classification%inner_wall_nodes)
        write(*,'("HC boundary faces: cyl_out=",I0," cyl_in=",I0," z-=",I0," z+=",I0," innerwall=",I0)') &
            classification%cyl_outer_faces%count, classification%cyl_inner_faces%count, &
            classification%axial_min_faces%count, classification%axial_max_faces%count, &
            classification%all_boundary_faces%count
    end subroutine report_hollowcylinder_summary

    subroutine report_box_summary(classification)
        type(BoxBoundaryClassification), intent(in) :: classification
        integer :: count_xmin, count_xmax, count_ymin, count_ymax, count_zmin, count_zmax

        count_xmin = size(classification%x_min_nodes)
        count_xmax = size(classification%x_max_nodes)
        count_ymin = size(classification%y_min_nodes)
        count_ymax = size(classification%y_max_nodes)
        count_zmin = size(classification%z_min_nodes)
        count_zmax = size(classification%z_max_nodes)

        write(*,'("Box boundary tolerance=",ES12.5)') classification%tolerance
        write(*,'("Box boundary nodes: x-=",I0," x+=",I0," y-=",I0," y+=",I0," z-=",I0," z+=",I0," innerwall=",I0)') &
            count_xmin, count_xmax, count_ymin, count_ymax, count_zmin, count_zmax, &
            size(classification%inner_wall_nodes)
        write(*,'("Box boundary faces: x-=",I0," x+=",I0," y-=",I0," y+=",I0," z-=",I0," z+=",I0," innerwall=",I0)') &
            classification%x_min_faces%count, classification%x_max_faces%count, classification%y_min_faces%count, &
            classification%y_max_faces%count, classification%z_min_faces%count, classification%z_max_faces%count, &
            classification%all_boundary_faces%count
    end subroutine report_box_summary

    subroutine write_hc_parametrizations(folder, classification, recorded_files, record_count)
        character(len=*), intent(in) :: folder
        type(HollowCylinderBoundaryClassification), intent(in) :: classification
        character(len=512), allocatable, intent(inout) :: recorded_files(:)
        integer, intent(inout) :: record_count

        call write_par_file(folder, "cyl_out.par", classification%cyl_outer_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "cyl_in.par", classification%cyl_inner_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "z-.par", classification%axial_min_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "z+.par", classification%axial_max_nodes, "Outflow", recorded_files, record_count)
        call write_par_file(folder, "innerwall.par", classification%inner_wall_nodes, "Wall", recorded_files, record_count)
    end subroutine write_hc_parametrizations

    subroutine write_box_parametrizations(folder, classification, recorded_files, record_count)
        character(len=*), intent(in) :: folder
        type(BoxBoundaryClassification), intent(in) :: classification
        character(len=512), allocatable, intent(inout) :: recorded_files(:)
        integer, intent(inout) :: record_count

        call write_par_file(folder, "x-.par", classification%x_min_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "x+.par", classification%x_max_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "y-.par", classification%y_min_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "y+.par", classification%y_max_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "z-.par", classification%z_min_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "z+.par", classification%z_max_nodes, "Outflow", recorded_files, record_count)
        call write_par_file(folder, "innerwall.par", classification%inner_wall_nodes, "Wall", recorded_files, record_count)
    end subroutine write_box_parametrizations

    subroutine write_par_file(folder, filename, nodes, keyword, recorded_files, record_count)
        character(len=*), intent(in) :: folder, filename
        integer, allocatable, intent(in) :: nodes(:)
        character(len=*), intent(in), optional :: keyword
        character(len=512), allocatable, intent(inout), optional :: recorded_files(:)
        integer, intent(inout), optional :: record_count
        integer :: unit, count, i
        character(len=512) :: filepath
        character(len=32) :: label

        call ensure_directory_exists(folder)
        call build_config_path(folder, filename, filepath)
        open(newunit=unit, file=trim(filepath), status="replace", action="write")
        if (allocated(nodes)) then
            count = size(nodes)
        else
            count = 0
        end if
        label = "Wall"
        if (present(keyword)) label = trim(keyword)
        write(unit,'(I0,1X,A)') count, trim(label)
        write(unit,'(A)') '" "'
        if (count > 0) then
            do i = 1, count
                write(unit,'(I0)') nodes(i)
            end do
        end if
        close(unit)
        if (present(recorded_files) .and. present(record_count)) then
        call append_file_record(recorded_files, record_count, trim(filename))
        end if
    end subroutine write_par_file

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
