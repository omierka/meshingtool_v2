program fortran_cgal_demo
    use iso_c_binding, only: c_ptr, c_null_char, c_associated, c_size_t, c_double, c_int, c_char
    use iso_fortran_env, only: output_unit, int64
    use mpi
    use bc_treatment, only: recompute_knpr_from_connectivity, HollowCylinderBoundaryClassification, &
        classify_hollowcylinder_boundaries, BoxBoundaryClassification, classify_box_boundaries, FaceList, &
        InflowBoundaryGroup
    use setupe3dfile_reader, only: MeshConfig, initialize_mesh_config, load_mesh_config, log_mesh_config, &
        ProcessParameters, initialize_process_parameters, load_process_parameters, log_process_inflows
    implicit none

    type :: WallFieldData
        character(len=64) :: name = ""
        integer(c_int), allocatable :: values(:)
    end type WallFieldData

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
    integer(c_int), allocatable :: hex_monitor(:)
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
    integer(c_int), allocatable :: filtered_monitor(:)
    integer(c_int), allocatable :: node_inflow_ids(:)
    type(WallFieldData), allocatable :: wall_fields(:)
    integer(c_int), allocatable :: filtered_kadj(:, :)
    real(c_double), allocatable :: filtered_volumes(:)
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
    type(ProcessParameters) :: process_params
    type(HollowCylinderBoundaryClassification) :: hc_boundary
    type(BoxBoundaryClassification) :: box_boundary
    type(FaceList) :: boundary_faces
    logical :: hc_summary_ready, box_summary_ready
    character(len=512) :: mesh_output_folder
    character(len=512) :: monitor_summary_file
    character(len=512) :: monitor_volume_file
    character(len=512), allocatable :: meshdir_files(:)
    integer :: meshdir_file_count
    integer :: mpi_rank, mpi_size, mpi_err
    integer :: worker_count, worker_rank
    logical :: is_master, is_worker
    integer, parameter :: PROGRESS_TAG = 9001
    integer, parameter :: MONITOR_BUCKETS = 4
    logical :: monitor_data_available, use_vtu_input
    integer :: monitor_hist_before(MONITOR_BUCKETS), monitor_hist_after(MONITOR_BUCKETS)
    real(c_double) :: monitor_volume_hist(MONITOR_BUCKETS)

    surface_file = "sphere.off"
    hex_file = "single.tri"
    hex_scale = 1.0_c_double
    output_folder= "."
    config_file = ""
    call initialize_mesh_config(mesh_config)
    call initialize_process_parameters(process_params)
    hc_summary_ready = .false.
    box_summary_ready = .false.
    meshdir_file_count = 0
    allocate(meshdir_files(0))
    monitor_data_available = .false.
    monitor_hist_before = 0
    monitor_hist_after = 0
    monitor_volume_hist = 0.0_c_double

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
        call load_process_parameters(trim(config_file), process_params)
        if (is_master) call log_process_inflows(process_params)
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

    use_vtu_input = file_has_extension(trim(hex_file), ".vtu")
    if (use_vtu_input) then
        call read_vtu_hex(trim(hex_file), hex_nel, hex_nvt, hex_nve, hex_nbct, hex_nee, hex_nae, header_line1, header_line2, &
                          hex_coords, hex_kvert, hex_knpr, hex_monitor, monitor_data_available)
        hex_nbct = 1
        hex_nve = 8
        hex_nee = 12
        hex_nae = 6
    else
        call read_single_tri(trim(hex_file), hex_nel, hex_nvt, hex_nve, hex_nbct, hex_nee, hex_nae, header_line1, header_line2, &
                             hex_coords, hex_kvert, hex_knpr)
        allocate(hex_monitor(hex_nel))
        hex_monitor = 0
        monitor_data_available = .false.
    end if
    hex_coords = hex_scale * hex_coords
    if (monitor_data_available .and. is_master) then
        call compute_monitor_histogram(hex_monitor, monitor_hist_before)
        call print_monitor_distribution("[MONITOR_BEFORE]", monitor_hist_before, hex_nel)
    end if
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
        if (monitor_data_available) then
            call compute_monitor_histogram_masked(hex_monitor, intersect_mask, monitor_hist_after)
            call print_monitor_distribution("[MONITOR_AFTER]", monitor_hist_after, intersect_count)
        end if
    end if

    if (is_master .and. intersect_count > 0) then
        call ensure_directory_exists(mesh_output_folder)
        call build_config_path(trim(mesh_output_folder), "Filtered.tri", filtered_hex_file)
        write(filtered_vtu_file(1:),'(A)') adjustl(trim(output_folder))//"/Filtered.vtu"
        write(filtered_tri_file(1:),'(A)') adjustl(trim(output_folder))//"/Filtered.tri"
        if (monitor_data_available) then
            call reduce_hex_mesh(intersect_mask, hex_coords, hex_kvert, hex_knpr, hex_nve, filtered_coords, &
                                 filtered_kvert, filtered_knpr, filtered_nel, filtered_nvt, monitor_in=hex_monitor, &
                                 out_monitor=filtered_monitor)
        else
            call reduce_hex_mesh(intersect_mask, hex_coords, hex_kvert, hex_knpr, hex_nve, filtered_coords, &
                                 filtered_kvert, filtered_knpr, filtered_nel, filtered_nvt)
        end if
        if (monitor_data_available) then
            call compute_hexahedron_volumes(filtered_coords, filtered_kvert, filtered_volumes)
            call compute_monitor_volume_histogram(filtered_monitor, filtered_volumes, monitor_volume_hist)
            if (allocated(filtered_volumes)) deallocate(filtered_volumes)
        end if
        if (allocated(filtered_kadj)) deallocate(filtered_kadj)
        allocate(filtered_kadj(6, filtered_nel))
        filtered_kadj = 0_c_int
        ! Build face-to-face adjacency before boundary reclassification for consistency checks.
        call build_hex_adjacency(filtered_kvert, filtered_kadj)
        call validate_hex_connectivity(filtered_kadj)
        call recompute_knpr_from_connectivity(filtered_kvert, filtered_knpr, boundary_faces)
        if (allocated(node_inflow_ids)) deallocate(node_inflow_ids)
        allocate(node_inflow_ids(filtered_nvt))
        node_inflow_ids = 0_c_int
        if (mesh_config%loaded) then
            call clear_wall_fields(wall_fields)
            select case (trim(mesh_config%mesh_type))
            case ("HollowCylinder","FullCylinder")
                call classify_hollowcylinder_boundaries(mesh_config, process_params, filtered_coords, filtered_kvert, &
                    filtered_knpr, boundary_faces, hc_boundary)
                hc_summary_ready = .true.
                call write_hc_parametrizations(mesh_output_folder, mesh_config, hc_boundary, meshdir_files, &
                    meshdir_file_count, filtered_nvt, wall_fields)
                call populate_hc_inflow_field(hc_boundary, node_inflow_ids)
            case ("Box")
                call classify_box_boundaries(mesh_config, process_params, filtered_coords, filtered_kvert, filtered_knpr, &
                    boundary_faces, &
                    box_boundary)
                box_summary_ready = .true.
                call write_box_parametrizations(mesh_output_folder, mesh_config, box_boundary, meshdir_files, &
                    meshdir_file_count, filtered_nvt, wall_fields)
                call populate_box_inflow_field(box_boundary, node_inflow_ids)
            end select
        end if
        call write_single_tri(filtered_tri_file, header_line1, header_line2, filtered_nel, filtered_nvt, hex_nbct, &
                              hex_nve, hex_nee, hex_nae, 1.0_c_double * filtered_coords, filtered_kvert, filtered_knpr)
        call write_single_tri(filtered_hex_file, header_line1, header_line2, filtered_nel, filtered_nvt, hex_nbct, &
                              hex_nve, hex_nee, hex_nae, 0.1_c_double * filtered_coords, filtered_kvert, filtered_knpr)
        call append_file_record(meshdir_files, meshdir_file_count, "Filtered.tri")
        if (monitor_data_available) then
            if (allocated(wall_fields)) then
                call write_vtu(filtered_vtu_file, filtered_coords, filtered_kvert, filtered_knpr, inflow_ids=node_inflow_ids, &
                    monitor_values=filtered_monitor, wall_fields=wall_fields)
            else
                call write_vtu(filtered_vtu_file, filtered_coords, filtered_kvert, filtered_knpr, inflow_ids=node_inflow_ids, &
                    monitor_values=filtered_monitor)
            end if
        else
            if (allocated(wall_fields)) then
                call write_vtu(filtered_vtu_file, filtered_coords, filtered_kvert, filtered_knpr, inflow_ids=node_inflow_ids, &
                    wall_fields=wall_fields)
            else
                call write_vtu(filtered_vtu_file, filtered_coords, filtered_kvert, filtered_knpr, inflow_ids=node_inflow_ids)
            end if
        end if
        call clear_wall_fields(wall_fields)
        call write_project_file(mesh_output_folder, meshdir_files, meshdir_file_count)
    end if

    if (is_master) elapsed_time = MPI_Wtime() - start_time

    call cgal_free_mesh(mesh_handle)

    if (is_master .and. monitor_data_available) then
        call ensure_directory_exists(output_folder)
        call build_config_path(trim(output_folder), "monitor_summary.txt", monitor_summary_file)
        call write_monitor_summary_file(trim(monitor_summary_file), monitor_hist_after, intersect_count)
        call build_config_path(trim(output_folder), "monitor_summary_volumetric.txt", monitor_volume_file)
        call write_monitor_volume_summary_file(trim(monitor_volume_file), monitor_volume_hist)
    end if

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
        integer :: inflow_count, inflow_idx
        character(len=32) :: inflow_label

        count_outer = size(classification%cyl_outer_nodes)
        count_inner = size(classification%cyl_inner_nodes)
        count_zmin = size(classification%axial_min_nodes)
        count_zmax = size(classification%axial_max_nodes)
        inflow_count = 0
        if (allocated(classification%inflow_groups)) inflow_count = size(classification%inflow_groups)

        write(*,'("HC boundary tolerance=",ES12.5)') classification%tolerance
        write(*,'("HC boundary nodes: cyl_out=",I0," cyl_in=",I0," z-=",I0," z+=",I0," innerwall=",I0)') &
            count_outer, count_inner, count_zmin, count_zmax, size(classification%inner_wall_nodes)
        write(*,'("HC boundary faces: cyl_out=",I0," cyl_in=",I0," z-=",I0," z+=",I0," innerwall=",I0)') &
            classification%cyl_outer_faces%count, classification%cyl_inner_faces%count, &
            classification%axial_min_faces%count, classification%axial_max_faces%count, &
            classification%all_boundary_faces%count
        if (inflow_count > 0) then
            write(*,'("HC inflow groups=",I0)') inflow_count
            do inflow_idx = 1, inflow_count
                if (len_trim(classification%inflow_groups(inflow_idx)%label) > 0) then
                    inflow_label = trim(classification%inflow_groups(inflow_idx)%label)
                else
                    inflow_label = "-"
                end if
                write(*,'("  inflow#",I0," type=",I0," faces=",I0," nodes=",I0," label=",A)') &
                    classification%inflow_groups(inflow_idx)%inflow_index, &
                    classification%inflow_groups(inflow_idx)%inflow_type, &
                    classification%inflow_groups(inflow_idx)%faces%count, &
                    size(classification%inflow_groups(inflow_idx)%nodes), trim(inflow_label)
            end do
        end if
    end subroutine report_hollowcylinder_summary

    subroutine report_box_summary(classification)
        type(BoxBoundaryClassification), intent(in) :: classification
        integer :: count_xmin, count_xmax, count_ymin, count_ymax, count_zmin, count_zmax
        integer :: inflow_count, inflow_idx
        character(len=32) :: inflow_label

        count_xmin = size(classification%x_min_nodes)
        count_xmax = size(classification%x_max_nodes)
        count_ymin = size(classification%y_min_nodes)
        count_ymax = size(classification%y_max_nodes)
        count_zmin = size(classification%z_min_nodes)
        count_zmax = size(classification%z_max_nodes)
        inflow_count = 0
        if (allocated(classification%inflow_groups)) inflow_count = size(classification%inflow_groups)

        write(*,'("Box boundary tolerance=",ES12.5)') classification%tolerance
        write(*,'("Box boundary nodes: x-=",I0," x+=",I0," y-=",I0," y+=",I0," z-=",I0," z+=",I0," innerwall=",I0)') &
            count_xmin, count_xmax, count_ymin, count_ymax, count_zmin, count_zmax, &
            size(classification%inner_wall_nodes)
        write(*,'("Box boundary faces: x-=",I0," x+=",I0," y-=",I0," y+=",I0," z-=",I0," z+=",I0," innerwall=",I0)') &
            classification%x_min_faces%count, classification%x_max_faces%count, classification%y_min_faces%count, &
            classification%y_max_faces%count, classification%z_min_faces%count, classification%z_max_faces%count, &
            classification%all_boundary_faces%count
        if (inflow_count > 0) then
            write(*,'("Box inflow groups=",I0)') inflow_count
            do inflow_idx = 1, inflow_count
                if (len_trim(classification%inflow_groups(inflow_idx)%label) > 0) then
                    inflow_label = trim(classification%inflow_groups(inflow_idx)%label)
                else
                    inflow_label = "-"
                end if
                write(*,'("  inflow#",I0," type=",I0," faces=",I0," nodes=",I0," label=",A)') &
                    classification%inflow_groups(inflow_idx)%inflow_index, &
                    classification%inflow_groups(inflow_idx)%inflow_type, &
                    classification%inflow_groups(inflow_idx)%faces%count, &
                    size(classification%inflow_groups(inflow_idx)%nodes), trim(inflow_label)
            end do
        end if
    end subroutine report_box_summary

    subroutine write_hc_parametrizations(folder, mesh_config, classification, recorded_files, record_count, &
        total_nodes, wall_fields)
        character(len=*), intent(in) :: folder
        type(MeshConfig), intent(in) :: mesh_config
        type(HollowCylinderBoundaryClassification), intent(in) :: classification
        character(len=512), allocatable, intent(inout) :: recorded_files(:)
        integer, intent(inout) :: record_count
        integer, intent(in) :: total_nodes
        type(WallFieldData), allocatable, intent(inout) :: wall_fields(:)
        integer :: inflow_idx
        character(len=64) :: inflow_filename
        character(len=32) :: inflow_label
        real(c_double) :: z_max_offset
        logical :: have_z_offset
        logical :: is_full_cyl

        is_full_cyl = (trim(mesh_config%mesh_type) == "FullCylinder")
        call write_par_file(folder, "cyl_out.par", classification%cyl_outer_nodes, "Wall", recorded_files, record_count)
        if (allocated(classification%cyl_outer_nodes)) then
            if (size(classification%cyl_outer_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_cyl_out", classification%cyl_outer_nodes)
        end if
        if (.not. is_full_cyl) then
            call write_par_file(folder, "cyl_in.par", classification%cyl_inner_nodes, "Wall", recorded_files, record_count)
            if (allocated(classification%cyl_inner_nodes)) then
                if (size(classification%cyl_inner_nodes) > 0) &
                    call append_wall_field(wall_fields, total_nodes, "Wall_cyl_in", classification%cyl_inner_nodes)
            end if
        end if
        call write_par_file(folder, "z-.par", classification%axial_min_nodes, "Wall", recorded_files, record_count)
        if (allocated(classification%axial_min_nodes)) then
            if (size(classification%axial_min_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_z-", classification%axial_min_nodes)
        end if
        have_z_offset = .false.
        if (mesh_config%cylinder%has_axial_start .and. mesh_config%cylinder%has_barrel_length) then
            z_max_offset = -1.0_c_double * (mesh_config%cylinder%axial_start + mesh_config%cylinder%barrel_length) * 0.1_c_double
            have_z_offset = .true.
        end if
        if (have_z_offset) then
            call write_par_file(folder, "z+.par", classification%axial_max_nodes, "Outflow", recorded_files, record_count, &
                offset_value=z_max_offset)
        else
            call write_par_file(folder, "z+.par", classification%axial_max_nodes, "Outflow", recorded_files, record_count)
        end if
        if (allocated(classification%axial_max_nodes)) then
            if (size(classification%axial_max_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_z+", classification%axial_max_nodes)
        end if
        call write_par_file(folder, "innerwall.par", classification%inner_wall_nodes, "Wall", recorded_files, record_count)
        if (allocated(classification%inner_wall_nodes)) then
            if (size(classification%inner_wall_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_innerwall", classification%inner_wall_nodes)
        end if
        if (allocated(classification%inflow_groups)) then
            do inflow_idx = 1, size(classification%inflow_groups)
                if (is_full_cyl) then
                    write(inflow_filename,'("fc_inflow_",I0,".par")') classification%inflow_groups(inflow_idx)%inflow_index
                else
                    write(inflow_filename,'("hc_inflow_",I0,".par")') classification%inflow_groups(inflow_idx)%inflow_index
                end if
                write(inflow_label,'("Inflow-",I0)') classification%inflow_groups(inflow_idx)%inflow_index
                call write_par_file(folder, trim(inflow_filename), classification%inflow_groups(inflow_idx)%nodes, &
                    trim(inflow_label), recorded_files, record_count)
            end do
        end if
    end subroutine write_hc_parametrizations

    subroutine write_box_parametrizations(folder, mesh_config, classification, recorded_files, record_count, &
        total_nodes, wall_fields)
        character(len=*), intent(in) :: folder
        type(MeshConfig), intent(in) :: mesh_config
        type(BoxBoundaryClassification), intent(in) :: classification
        character(len=512), allocatable, intent(inout) :: recorded_files(:)
        integer, intent(inout) :: record_count
        integer, intent(in) :: total_nodes
        type(WallFieldData), allocatable, intent(inout) :: wall_fields(:)
        integer :: inflow_idx
        character(len=64) :: inflow_filename
        character(len=32) :: inflow_label
        real(c_double) :: z_max_offset
        logical :: have_z_offset

        call write_par_file(folder, "x-.par", classification%x_min_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "x+.par", classification%x_max_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "y-.par", classification%y_min_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "y+.par", classification%y_max_nodes, "Wall", recorded_files, record_count)
        call write_par_file(folder, "z-.par", classification%z_min_nodes, "Wall", recorded_files, record_count)
        if (allocated(classification%x_min_nodes)) then
            if (size(classification%x_min_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_x-", classification%x_min_nodes)
        end if
        if (allocated(classification%x_max_nodes)) then
            if (size(classification%x_max_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_x+", classification%x_max_nodes)
        end if
        if (allocated(classification%y_min_nodes)) then
            if (size(classification%y_min_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_y-", classification%y_min_nodes)
        end if
        if (allocated(classification%y_max_nodes)) then
            if (size(classification%y_max_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_y+", classification%y_max_nodes)
        end if
        if (allocated(classification%z_min_nodes)) then
            if (size(classification%z_min_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_z-", classification%z_min_nodes)
        end if
        have_z_offset = .false.
        if (mesh_config%box%has_geometry_start .and. mesh_config%box%has_geometry_length) then
            z_max_offset = -1.0_c_double * &
                (mesh_config%box%geometry_start(3) + mesh_config%box%geometry_length(3)) * 0.1_c_double
            have_z_offset = .true.
        end if
        if (have_z_offset) then
            call write_par_file(folder, "z+.par", classification%z_max_nodes, "Outflow", recorded_files, record_count, &
                offset_value=z_max_offset)
        else
            call write_par_file(folder, "z+.par", classification%z_max_nodes, "Outflow", recorded_files, record_count)
        end if
        call write_par_file(folder, "innerwall.par", classification%inner_wall_nodes, "Wall", recorded_files, record_count)
        if (allocated(classification%z_max_nodes)) then
            if (size(classification%z_max_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_z+", classification%z_max_nodes)
        end if
        if (allocated(classification%inner_wall_nodes)) then
            if (size(classification%inner_wall_nodes) > 0) &
                call append_wall_field(wall_fields, total_nodes, "Wall_innerwall", classification%inner_wall_nodes)
        end if
        if (allocated(classification%inflow_groups)) then
            do inflow_idx = 1, size(classification%inflow_groups)
                write(inflow_filename,'("inflow_",I0,".par")') classification%inflow_groups(inflow_idx)%inflow_index
                write(inflow_label,'("Inflow-",I0)') classification%inflow_groups(inflow_idx)%inflow_index
                call write_par_file(folder, trim(inflow_filename), classification%inflow_groups(inflow_idx)%nodes, &
                    trim(inflow_label), recorded_files, record_count)
            end do
        end if
    end subroutine write_box_parametrizations

    subroutine populate_hc_inflow_field(classification, inflow_ids)
        type(HollowCylinderBoundaryClassification), intent(in) :: classification
        integer(c_int), intent(inout) :: inflow_ids(:)

        if (allocated(classification%inflow_groups)) then
            call assign_inflow_groups_to_nodes(classification%inflow_groups, inflow_ids)
        end if
        if (allocated(classification%axial_max_nodes)) then
            call stamp_nodes_with_value(classification%axial_max_nodes, -1_c_int, inflow_ids)
        end if
    end subroutine populate_hc_inflow_field

    subroutine populate_box_inflow_field(classification, inflow_ids)
        type(BoxBoundaryClassification), intent(in) :: classification
        integer(c_int), intent(inout) :: inflow_ids(:)

        if (allocated(classification%inflow_groups)) then
            call assign_inflow_groups_to_nodes(classification%inflow_groups, inflow_ids)
        end if
        if (allocated(classification%z_max_nodes)) then
            call stamp_nodes_with_value(classification%z_max_nodes, -1_c_int, inflow_ids)
        end if
    end subroutine populate_box_inflow_field

    subroutine assign_inflow_groups_to_nodes(groups, inflow_ids)
        type(InflowBoundaryGroup), intent(in) :: groups(:)
        integer(c_int), intent(inout) :: inflow_ids(:)
        integer :: inflow_idx
        integer(c_int) :: inflow_value

        if (size(groups) == 0) return
        do inflow_idx = 1, size(groups)
            if (.not. allocated(groups(inflow_idx)%nodes)) cycle
            inflow_value = int(groups(inflow_idx)%inflow_index, kind=c_int)
            call stamp_nodes_with_value(groups(inflow_idx)%nodes, inflow_value, inflow_ids)
        end do
    end subroutine assign_inflow_groups_to_nodes

    subroutine stamp_nodes_with_value(nodes, value, inflow_ids)
        integer, intent(in) :: nodes(:)
        integer(c_int), intent(in) :: value
        integer(c_int), intent(inout) :: inflow_ids(:)
        integer :: idx, node_id, limit

        limit = size(inflow_ids)
        do idx = 1, size(nodes)
            node_id = nodes(idx)
            if (node_id < 1 .or. node_id > limit) cycle
            inflow_ids(node_id) = value
        end do
    end subroutine stamp_nodes_with_value

    subroutine write_par_file(folder, filename, nodes, keyword, recorded_files, record_count, offset_value)
        character(len=*), intent(in) :: folder, filename
        integer, allocatable, intent(in) :: nodes(:)
        character(len=*), intent(in), optional :: keyword
        character(len=512), allocatable, intent(inout), optional :: recorded_files(:)
        integer, intent(inout), optional :: record_count
        real(c_double), intent(in), optional :: offset_value
        integer :: unit, count, i
        character(len=512) :: filepath
        character(len=32) :: label
        character(len=64) :: header_line

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
        if (present(offset_value)) then
            write(header_line,'("4 0.0 0.0 1.0 ",F16.6)') offset_value
            write(unit,'(A)') '"'//trim(adjustl(header_line))//'"'
        else
            write(unit,'(A)') '"1 0.0"'
        end if
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

    subroutine append_wall_field(fields, total_nodes, field_label, nodes)
        type(WallFieldData), allocatable, intent(inout) :: fields(:)
        integer, intent(in) :: total_nodes
        character(len=*), intent(in) :: field_label
        integer, intent(in) :: nodes(:)
        type(WallFieldData), allocatable :: temp(:)
        integer :: old_count, idx, node_idx, node_id

        if (total_nodes <= 0) return
        if (size(nodes) <= 0) return

        old_count = 0
        if (allocated(fields)) then
            old_count = size(fields)
            if (old_count > 0) then
                allocate(temp(old_count))
                temp = fields
            end if
            call clear_wall_fields(fields)
            allocate(fields(old_count + 1))
            if (old_count > 0) then
                fields(1:old_count) = temp
                deallocate(temp)
            end if
        else
            allocate(fields(1))
        end if
        idx = old_count + 1
        fields(idx)%name = trim(field_label)
        if (allocated(fields(idx)%values)) deallocate(fields(idx)%values)
        allocate(fields(idx)%values(total_nodes))
        fields(idx)%values = 0_c_int
        do node_idx = 1, size(nodes)
            node_id = nodes(node_idx)
            if (node_id >= 1 .and. node_id <= total_nodes) fields(idx)%values(node_id) = 1_c_int
        end do
    end subroutine append_wall_field

    subroutine clear_wall_fields(fields)
        type(WallFieldData), allocatable, intent(inout) :: fields(:)
        integer :: idx
        if (.not.allocated(fields)) return
        do idx = 1, size(fields)
            if (allocated(fields(idx)%values)) deallocate(fields(idx)%values)
            fields(idx)%name = ""
        end do
        deallocate(fields)
    end subroutine clear_wall_fields

    subroutine compute_monitor_histogram(monitor_values, histogram)
        integer(c_int), intent(in) :: monitor_values(:)
        integer, intent(out) :: histogram(MONITOR_BUCKETS)
        integer :: idx, bucket

        histogram = 0
        do idx = 1, size(monitor_values)
            bucket = int(monitor_values(idx))
            if (bucket >= 0 .and. bucket < MONITOR_BUCKETS) then
                histogram(bucket + 1) = histogram(bucket + 1) + 1
            end if
        end do
    end subroutine compute_monitor_histogram

    subroutine compute_monitor_histogram_masked(monitor_values, mask, histogram)
        integer(c_int), intent(in) :: monitor_values(:)
        logical, intent(in) :: mask(:)
        integer, intent(out) :: histogram(MONITOR_BUCKETS)
        integer :: idx, bucket

        if (size(monitor_values) /= size(mask)) error stop "Monitor histogram mask mismatch."
        histogram = 0
        do idx = 1, size(monitor_values)
            if (.not. mask(idx)) cycle
            bucket = int(monitor_values(idx))
            if (bucket >= 0 .and. bucket < MONITOR_BUCKETS) then
                histogram(bucket + 1) = histogram(bucket + 1) + 1
            end if
        end do
    end subroutine compute_monitor_histogram_masked

    subroutine compute_monitor_volume_histogram(monitor_values, volumes, histogram)
        integer(c_int), intent(in) :: monitor_values(:)
        real(c_double), intent(in) :: volumes(:)
        real(c_double), intent(out) :: histogram(MONITOR_BUCKETS)
        integer :: idx, bucket

        if (size(monitor_values) /= size(volumes)) error stop "Monitor volume histogram size mismatch."
        histogram = 0.0_c_double
        do idx = 1, size(monitor_values)
            bucket = int(monitor_values(idx))
            if (bucket >= 0 .and. bucket < MONITOR_BUCKETS) then
                histogram(bucket + 1) = histogram(bucket + 1) + max(volumes(idx), 0.0_c_double)
            end if
        end do
    end subroutine compute_monitor_volume_histogram

    subroutine print_monitor_distribution(label, histogram, total_count)
        character(len=*), intent(in) :: label
        integer, intent(in) :: histogram(MONITOR_BUCKETS)
        integer, intent(in) :: total_count
        real(c_double) :: pct(MONITOR_BUCKETS)
        real(c_double) :: denom
        integer :: bucket

        denom = real(max(1, total_count), kind=c_double)
        do bucket = 1, MONITOR_BUCKETS
            pct(bucket) = 100.0_c_double * real(histogram(bucket), kind=c_double) / denom
        end do
        write (*,'(A,1X,"0=[",F6.2,"%]",1X,"1=[",F6.2,"%]",1X,"2=[",F6.2,"%]",1X,"3=[",F6.2,"%]")') &
            trim(label), pct(1), pct(2), pct(3), pct(4)
    end subroutine print_monitor_distribution

    subroutine write_monitor_summary_file(filepath, histogram, total_count)
        character(len=*), intent(in) :: filepath
        integer, intent(in) :: histogram(MONITOR_BUCKETS)
        integer, intent(in) :: total_count
        integer :: unit
        real(c_double) :: pct(MONITOR_BUCKETS), denom

        denom = real(max(1, total_count), kind=c_double)
        pct = 0.0_c_double
        pct = 100.0_c_double * real(histogram, kind=c_double) / denom
        open(newunit=unit, file=trim(filepath), status="replace", action="write")
        write(unit,'("MONITOR_AFTER_COUNT ", "0=[",F5.1,"%]",1X,"1=[",F5.1,"%]",1X,"2=[",F5.1,"%]",1X,"3=[",F5.1,"%]")') &
            pct(1), pct(2), pct(3), pct(4)
        close(unit)
    end subroutine write_monitor_summary_file

    subroutine write_monitor_volume_summary_file(filepath, volume_histogram)
        character(len=*), intent(in) :: filepath
        real(c_double), intent(in) :: volume_histogram(MONITOR_BUCKETS)
        integer :: unit
        real(c_double) :: vol_pct(MONITOR_BUCKETS), volume_total

        volume_total = sum(volume_histogram)
        if (volume_total <= 0.0_c_double) then
            vol_pct = 0.0_c_double
        else
            vol_pct = 100.0_c_double * volume_histogram / volume_total
        end if
        open(newunit=unit, file=trim(filepath), status="replace", action="write")
        write(unit,'("MONITOR_AFTER_VOLUME ", "0=[",F5.1,"%]",1X,"1=[",F5.1,"%]",1X,"2=[",F5.1,"%]",1X,"3=[",F5.1,"%]")') &
            vol_pct(1), vol_pct(2), vol_pct(3), vol_pct(4)
        close(unit)
    end subroutine write_monitor_volume_summary_file

    pure function to_lower_char(ch) result(lower)
        character(len=1), intent(in) :: ch
        character(len=1) :: lower
        integer :: code

        lower = ch
        code = iachar(ch)
        if (code >= iachar('A') .and. code <= iachar('Z')) lower = achar(code + 32)
    end function to_lower_char

    logical function file_has_extension(path, extension) result(has_ext)
        character(len=*), intent(in) :: path, extension
        integer :: path_len, ext_len, idx, path_pos

        has_ext = .false.
        path_len = len_trim(path)
        ext_len = len_trim(extension)
        if (ext_len == 0 .or. path_len < ext_len) return
        path_pos = path_len - ext_len + 1
        has_ext = .true.
        do idx = 1, ext_len
            if (to_lower_char(path(path_pos:path_pos)) /= to_lower_char(extension(idx:idx))) then
                has_ext = .false.
                exit
            end if
            path_pos = path_pos + 1
        end do
    end function file_has_extension

    integer function extract_attribute_int(line, attribute) result(value)
        character(len=*), intent(in) :: line, attribute
        integer :: attr_pos, closing_pos, ios
        character(len=256) :: buffer

        value = 0
        buffer = ""
        attr_pos = index(line, trim(attribute))
        if (attr_pos <= 0) error stop "Missing VTU attribute."
        attr_pos = attr_pos + len_trim(attribute)
        closing_pos = index(line(attr_pos:), '"')
        if (closing_pos <= 0) error stop "Malformed VTU attribute."
        buffer = line(attr_pos:attr_pos + closing_pos - 2)
        read(buffer, *, iostat=ios) value
        if (ios /= 0) error stop "Failed to parse VTU attribute."
    end function extract_attribute_int

    subroutine read_vtu_hex(filename, nel, nvt, nve, nbct, nee, nae, header1, header2, coords, kvert, knpr, monitor, &
        monitor_available)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: nel, nvt, nve, nbct, nee, nae
        character(len=*), intent(out) :: header1, header2
        real(c_double), allocatable, intent(out) :: coords(:, :)
        integer(c_int), allocatable, intent(out) :: kvert(:, :)
        integer(c_int), allocatable, intent(out) :: knpr(:)
        integer(c_int), allocatable, intent(out) :: monitor(:)
        logical, intent(out) :: monitor_available

        integer :: unit, ios, cell_idx, local_idx, flat_idx
        character(len=1024) :: line
        logical :: in_pointdata, in_celldata, in_points, in_cells
        logical :: counts_ready, coords_loaded, knpr_loaded, monitor_loaded, conn_loaded
        real(c_double), allocatable :: coord_flat(:)
        integer(c_int), allocatable :: connectivity(:), offsets(:), types(:)

        nel = 0
        nvt = 0
        nve = 8
        nbct = 0
        nee = 0
        nae = 0
        header1 = "VTU INPUT"
        header2 = "Converted from VTU"
        monitor_available = .false.
        in_pointdata = .false.
        in_celldata = .false.
        in_points = .false.
        in_cells = .false.
        counts_ready = .false.
        coords_loaded = .false.
        knpr_loaded = .false.
        monitor_loaded = .false.
        conn_loaded = .false.

        open(newunit=unit, file=trim(filename), status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "Failed to open VTU mesh."

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, "<Piece") > 0) then
                nvt = extract_attribute_int(line, 'NumberOfPoints="')
                nel = extract_attribute_int(line, 'NumberOfCells="')
                if (nvt <= 0 .or. nel <= 0) error stop "Invalid VTU mesh sizes."
                if (.not. counts_ready) then
                    allocate(coords(3, nvt))
                    allocate(kvert(nve, nel))
                    allocate(knpr(nvt))
                    allocate(monitor(nel))
                    monitor = 0
                    counts_ready = .true.
                end if
            else if (index(line, "<PointData") > 0) then
                in_pointdata = .true.
            else if (index(line, "</PointData") > 0) then
                in_pointdata = .false.
            else if (index(line, "<CellData") > 0) then
                in_celldata = .true.
            else if (index(line, "</CellData") > 0) then
                in_celldata = .false.
            else if (index(line, "<Points") > 0) then
                in_points = .true.
            else if (index(line, "</Points") > 0) then
                in_points = .false.
            else if (index(line, "<Cells") > 0) then
                in_cells = .true.
            else if (index(line, "</Cells") > 0) then
                in_cells = .false.
            else if (in_pointdata .and. index(line, "<DataArray") > 0 .and. index(line, 'Name="KNPR"') > 0) then
                if (.not. counts_ready) error stop "VTU PointData encountered before Piece definition."
                call read_vtu_int_array(unit, nvt, knpr)
                knpr_loaded = .true.
            else if (in_celldata .and. index(line, "<DataArray") > 0 .and. index(line, 'Name="monitor"') > 0) then
                if (.not. counts_ready) error stop "VTU CellData encountered before Piece definition."
                call read_vtu_int_array(unit, nel, monitor)
                monitor_loaded = .true.
                monitor_available = .true.
            else if (in_points .and. index(line, "<DataArray") > 0 .and. index(line, 'NumberOfComponents="3"') > 0) then
                if (.not. counts_ready) error stop "VTU Points encountered before Piece definition."
                if (.not. allocated(coord_flat)) allocate(coord_flat(3 * nvt))
                call read_vtu_real_array(unit, 3 * nvt, coord_flat)
                coords_loaded = .true.
            else if (in_cells .and. index(line, "<DataArray") > 0 .and. index(line, 'Name="connectivity"') > 0) then
                if (.not. counts_ready) error stop "VTU Cells encountered before Piece definition."
                if (.not. allocated(connectivity)) allocate(connectivity(nve * nel))
                call read_vtu_int_array(unit, nve * nel, connectivity)
                conn_loaded = .true.
            else if (in_cells .and. index(line, "<DataArray") > 0 .and. index(line, 'Name="offsets"') > 0) then
                if (.not. counts_ready) error stop "VTU offsets encountered before Piece definition."
                if (.not. allocated(offsets)) allocate(offsets(nel))
                call read_vtu_int_array(unit, nel, offsets)
            else if (in_cells .and. index(line, "<DataArray") > 0 .and. index(line, 'Name="types"') > 0) then
                if (.not. counts_ready) error stop "VTU types encountered before Piece definition."
                if (.not. allocated(types)) allocate(types(nel))
                call read_vtu_int_array(unit, nel, types)
            end if
        end do

        close(unit)

        if (.not. counts_ready) error stop "VTU file missing Piece definition."
        if (.not. coords_loaded) error stop "VTU file missing coordinate data."
        if (.not. knpr_loaded) then
            knpr = 0
        end if
        if (.not. monitor_loaded) then
            monitor = 0
            monitor_available = .false.
        end if
        if (.not. conn_loaded) error stop "VTU file missing connectivity."

        if (allocated(types)) then
            do cell_idx = 1, nel
                if (int(types(cell_idx)) /= 12) then
                    error stop "Unsupported VTU cell type detected."
                end if
            end do
        end if

        flat_idx = 0
        do cell_idx = 1, nel
            do local_idx = 1, nve
                flat_idx = flat_idx + 1
                kvert(local_idx, cell_idx) = int(connectivity(flat_idx) + 1, kind=c_int)
            end do
        end do
        if (allocated(coord_flat)) then
            do cell_idx = 1, nvt
                coords(1, cell_idx) = coord_flat(3 * (cell_idx - 1) + 1)
                coords(2, cell_idx) = coord_flat(3 * (cell_idx - 1) + 2)
                coords(3, cell_idx) = coord_flat(3 * (cell_idx - 1) + 3)
            end do
            deallocate(coord_flat)
        end if
        if (allocated(connectivity)) deallocate(connectivity)
        if (allocated(offsets)) deallocate(offsets)
        if (allocated(types)) deallocate(types)
    end subroutine read_vtu_hex

    subroutine read_vtu_real_array(unit, count, values)
        integer, intent(in) :: unit, count
        real(c_double), intent(out) :: values(:)
        integer :: ios, idx
        character(len=1024) :: line

        if (count <= 0) then
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) error stop "Unexpected end of VTU data array."
            return
        end if
        if (size(values) < count) error stop "Insufficient VTU real buffer."
        read(unit, *)(values(idx), idx = 1, count)
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) error stop "Missing VTU DataArray terminator."
    end subroutine read_vtu_real_array

    subroutine read_vtu_int_array(unit, count, values)
        integer, intent(in) :: unit, count
        integer(c_int), intent(out) :: values(:)
        integer :: ios, idx
        character(len=1024) :: line

        if (count <= 0) then
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) error stop "Unexpected end of VTU data array."
            return
        end if
        if (size(values) < count) error stop "Insufficient VTU integer buffer."
        read(unit, *)(values(idx), idx = 1, count)
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) error stop "Missing VTU DataArray terminator."
    end subroutine read_vtu_int_array

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

    subroutine reduce_hex_mesh(intersect_mask, coords, kvert, knpr, nve, out_coords, out_kvert, out_knpr, out_nel, out_nvt, &
        monitor_in, out_monitor)
        logical, intent(in) :: intersect_mask(:)
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: kvert(:, :)
        integer(c_int), intent(in) :: knpr(:)
        integer, intent(in) :: nve
        real(c_double), allocatable, intent(out) :: out_coords(:, :)
        integer(c_int), allocatable, intent(out) :: out_kvert(:, :)
        integer(c_int), allocatable, intent(out) :: out_knpr(:)
        integer, intent(out) :: out_nel, out_nvt
        integer(c_int), intent(in), optional :: monitor_in(:)
        integer(c_int), allocatable, intent(out), optional :: out_monitor(:)

        integer :: hex_nel, hex_nvt
        logical, allocatable :: vertex_used(:)
        integer, allocatable :: vertex_remap(:)
        integer, allocatable :: temp_kvert(:, :)
        integer(c_int), allocatable :: temp_monitor(:)
        integer :: elem_idx, local_idx, new_elem_idx, vid, new_vid
        logical :: keep_monitor

        hex_nel = size(kvert, 2)
        hex_nvt = size(coords, 2)
        out_nel = count(intersect_mask)
        if (out_nel == 0) then
            out_nvt = 0
            allocate(out_coords(3, 0))
            allocate(out_kvert(nve, 0))
            allocate(out_knpr(0))
            if (present(out_monitor)) then
                if (.not. present(monitor_in)) error stop "Monitor output requested without input data."
                allocate(out_monitor(0))
            end if
            return
        end if

        allocate(vertex_used(hex_nvt))
        vertex_used = .false.
        allocate(temp_kvert(nve, out_nel))
        keep_monitor = present(monitor_in) .and. present(out_monitor)
        if (keep_monitor) then
            if (size(monitor_in) /= hex_nel) error stop "Monitor input size mismatch."
            allocate(temp_monitor(out_nel))
        end if

        new_elem_idx = 0
        do elem_idx = 1, hex_nel
            if (intersect_mask(elem_idx)) then
                new_elem_idx = new_elem_idx + 1
                do local_idx = 1, nve
                    vid = int(kvert(local_idx, elem_idx))
                    vertex_used(vid) = .true.
                    temp_kvert(local_idx, new_elem_idx) = vid
                end do
                if (keep_monitor) temp_monitor(new_elem_idx) = monitor_in(elem_idx)
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

        if (keep_monitor) then
            allocate(out_monitor(out_nel))
            out_monitor = temp_monitor
            deallocate(temp_monitor)
        end if

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

    subroutine write_vtu(filename, coords, kvert, knpr, inflow_ids, monitor_values, wall_fields)
        character(len=*), intent(in) :: filename
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: kvert(:, :)
        integer(c_int), intent(in) :: knpr(:)
        integer(c_int), intent(in), optional :: inflow_ids(:)
        integer(c_int), intent(in), optional :: monitor_values(:)
        type(WallFieldData), intent(in), optional :: wall_fields(:)

        integer :: unit, nvt, nel, nve
        integer :: i, j
        integer :: field_idx
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
        if (present(inflow_ids)) then
            if (size(inflow_ids) /= nvt) error stop "VTU inflow array size mismatch."
            write(unit, '(A)') '        <DataArray type="Int32" Name="InflowId" format="ascii">'
            write(unit, '(6(1X,I8))') (int(inflow_ids(i), kind=4), i = 1, nvt)
            write(unit, '(A)') '        </DataArray>'
        end if
        if (present(wall_fields)) then
            do field_idx = 1, size(wall_fields)
                if (.not.allocated(wall_fields(field_idx)%values)) cycle
                if (size(wall_fields(field_idx)%values) /= nvt) cycle
                if (len_trim(wall_fields(field_idx)%name) == 0) cycle
                write(unit, '(A)') '        <DataArray type="Int32" Name="'//trim(wall_fields(field_idx)%name)//'" format="ascii">'
                write(unit, '(6(1X,I8))') (int(wall_fields(field_idx)%values(i), kind=4), i = 1, nvt)
                write(unit, '(A)') '        </DataArray>'
            end do
        end if
        write(unit, '(A)') '      </PointData>'

        if (present(monitor_values)) then
            if (size(monitor_values) /= nel) error stop "VTU monitor data size mismatch."
            write(unit, '(A)') '      <CellData>'
            write(unit, '(A)') '        <DataArray type="Int32" Name="monitor" format="ascii">'
            write(unit, '(*(1X,I8))') (int(monitor_values(i), kind=4), i = 1, nel)
            write(unit, '(A)') '        </DataArray>'
            write(unit, '(A)') '      </CellData>'
        else
            write(unit, '(A)') '      <CellData/>'
        end if

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

    subroutine compute_hexahedron_volumes(coords, kvert, volumes)
        real(c_double), intent(in) :: coords(:, :)
        integer(c_int), intent(in) :: kvert(:, :)
        real(c_double), allocatable, intent(out) :: volumes(:)

        integer, parameter :: local_faces(4, 6) = reshape([ &
            1, 2, 3, 4, &
            5, 6, 7, 8, &
            1, 2, 6, 5, &
            2, 3, 7, 6, &
            3, 4, 8, 7, &
            4, 1, 5, 8], [4, 6])
        integer :: nel, nve, cell_idx, face_idx
        integer :: v1, v2, v3, v4, local_idx
        integer :: vertex_id
        real(c_double) :: cell_coords(3, 8)
        real(c_double) :: centroid(3)
        real(c_double) :: volume_sum

        nve = size(kvert, 1)
        if (nve /= 8) error stop "Hex volume computation expects 8-node elements."
        nel = size(kvert, 2)
        allocate(volumes(nel))
        volumes = 0.0_c_double
        if (nel == 0) return

        do cell_idx = 1, nel
            centroid = 0.0_c_double
            do local_idx = 1, nve
                vertex_id = int(kvert(local_idx, cell_idx))
                if (vertex_id < 1 .or. vertex_id > size(coords, 2)) error stop "Invalid vertex index in hex volume."
                cell_coords(:, local_idx) = coords(:, vertex_id)
                centroid = centroid + cell_coords(:, local_idx)
            end do
            centroid = centroid / 8.0_c_double
            volume_sum = 0.0_c_double
            do face_idx = 1, size(local_faces, 2)
                v1 = local_faces(1, face_idx)
                v2 = local_faces(2, face_idx)
                v3 = local_faces(3, face_idx)
                v4 = local_faces(4, face_idx)
                volume_sum = volume_sum + tetrahedron_volume(centroid, cell_coords(:, v1), cell_coords(:, v2), cell_coords(:, v3))
                volume_sum = volume_sum + tetrahedron_volume(centroid, cell_coords(:, v1), cell_coords(:, v3), cell_coords(:, v4))
            end do
            volumes(cell_idx) = volume_sum
        end do
    end subroutine compute_hexahedron_volumes

    pure function tetrahedron_volume(a, b, c, d) result(volume)
        real(c_double), intent(in) :: a(3), b(3), c(3), d(3)
        real(c_double) :: volume
        real(c_double) :: ab(3), ac(3), ad(3), cross_prod(3)

        ab = b - a
        ac = c - a
        ad = d - a
        cross_prod(1) = ac(2) * ad(3) - ac(3) * ad(2)
        cross_prod(2) = ac(3) * ad(1) - ac(1) * ad(3)
        cross_prod(3) = ac(1) * ad(2) - ac(2) * ad(1)
        volume = abs(ab(1) * cross_prod(1) + ab(2) * cross_prod(2) + ab(3) * cross_prod(3)) / 6.0_c_double
    end function tetrahedron_volume

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

    subroutine build_hex_adjacency(kvert, kadj)
        integer(c_int), intent(in) :: kvert(:, :)
        integer(c_int), intent(inout) :: kadj(:, :)

        integer, parameter :: faces_per_hex = 6
        integer, parameter :: face_pattern(4, faces_per_hex) = reshape([ &
            1, 2, 3, 4, &
            1, 2, 6, 5, &
            2, 3, 7, 6, &
            3, 4, 8, 7, &
            4, 1, 5, 8, &
            5, 6, 7, 8], [4, faces_per_hex])
        integer :: nel, nfaces, elem_idx, face_idx, entry_count
        integer :: bucket_idx, bucket_entry, owner_elem, owner_face
        integer :: hash_size, conflict_count
        integer(c_int) :: local_vertices(4), sorted_vertices(4)
        integer, allocatable :: hash_head(:), hash_next(:)
        integer, allocatable :: face_owner(:), face_local(:)
        integer(c_int), allocatable :: stored_vertices(:, :)
        logical :: match_found

        nel = size(kvert, 2)
        if (nel <= 0) return

        nfaces = faces_per_hex * nel
        if (nfaces <= 0) return

        hash_size = max(1, 2 * nfaces + 1)
        allocate(hash_head(hash_size))
        allocate(hash_next(nfaces))
        allocate(face_owner(nfaces))
        allocate(face_local(nfaces))
        allocate(stored_vertices(4, nfaces))
        hash_head = 0
        hash_next = 0
        face_owner = 0
        face_local = 0
        stored_vertices = 0_c_int
        entry_count = 0
        conflict_count = 0

        do elem_idx = 1, nel
            do face_idx = 1, faces_per_hex
                local_vertices = kvert(face_pattern(:, face_idx), elem_idx)
                sorted_vertices = local_vertices
                call sort_face_vertices(sorted_vertices)
                bucket_idx = face_hash(sorted_vertices, hash_size)
                bucket_entry = hash_head(bucket_idx)
                match_found = .false.
                do while (bucket_entry /= 0)
                    if (all(stored_vertices(:, bucket_entry) == sorted_vertices)) then
                        owner_elem = face_owner(bucket_entry)
                        owner_face = face_local(bucket_entry)
                        match_found = .true.
                        if (kadj(owner_face, owner_elem) /= 0_c_int .and. &
                            int(kadj(owner_face, owner_elem)) /= elem_idx) then
                            conflict_count = conflict_count + 1
                        else
                            kadj(face_idx, elem_idx) = int(owner_elem, kind=c_int)
                            kadj(owner_face, owner_elem) = int(elem_idx, kind=c_int)
                        end if
                        exit
                    end if
                    bucket_entry = hash_next(bucket_entry)
                end do
                if (.not. match_found) then
                    entry_count = entry_count + 1
                    if (entry_count > nfaces) then
                        error stop "Adjacency accumulator overflow."
                    end if
                    stored_vertices(:, entry_count) = sorted_vertices
                    face_owner(entry_count) = elem_idx
                    face_local(entry_count) = face_idx
                    hash_next(entry_count) = hash_head(bucket_idx)
                    hash_head(bucket_idx) = entry_count
                end if
            end do
        end do

        if (conflict_count > 0 .and. is_master) then
            write (*,'(A,I0)') '[ADJ_WARN] faces shared by more than 2 elements: ', conflict_count
        end if

        deallocate(hash_head, hash_next, face_owner, face_local, stored_vertices)
    end subroutine build_hex_adjacency

    subroutine validate_hex_connectivity(kadj)
        integer(c_int), intent(in) :: kadj(:, :)

        integer :: nel, faces_per_hex
        logical, allocatable :: visited(:)
        integer, allocatable :: queue(:)
        integer, allocatable :: component_sizes(:)
        integer :: elem_idx, face_idx, neighbor_elem
        integer :: head, tail, current_elem
        integer :: component_count, comp_idx

        nel = size(kadj, 2)
        faces_per_hex = size(kadj, 1)
        if (nel <= 0) return

        allocate(visited(nel))
        allocate(queue(nel))
        allocate(component_sizes(nel))
        visited = .false.
        component_sizes = 0
        component_count = 0

        do elem_idx = 1, nel
            if (.not. visited(elem_idx)) then
                component_count = component_count + 1
                head = 1
                tail = 1
                queue(tail) = elem_idx
                visited(elem_idx) = .true.
                do while (head <= tail)
                    current_elem = queue(head)
                    head = head + 1
                    component_sizes(component_count) = component_sizes(component_count) + 1
                    do face_idx = 1, faces_per_hex
                        neighbor_elem = int(kadj(face_idx, current_elem))
                        if (neighbor_elem > 0 .and. neighbor_elem <= nel) then
                            if (.not. visited(neighbor_elem)) then
                                tail = tail + 1
                                queue(tail) = neighbor_elem
                                visited(neighbor_elem) = .true.
                            end if
                        end if
                    end do
                end do
            end if
        end do

        if (is_master) then
            write (*,'(A)') '[CONNECTIVITY] Checking face-adjacent regions...'
            write (*,'(A,I0,A)') '[CONNECTIVITY] Identified ', component_count, ' independent subregions.'
            do comp_idx = 1, component_count
                write (*,'(A,I0,A,I0)') '[CONNECTIVITY] Subregion ', comp_idx, ': ', component_sizes(comp_idx), &
                    ' elements'
            end do
        end if

        deallocate(visited, queue, component_sizes)
    end subroutine validate_hex_connectivity

    subroutine sort_face_vertices(vertices)
        integer(c_int), intent(inout) :: vertices(4)
        integer(c_int) :: key
        integer :: i, j

        do i = 2, 4
            key = vertices(i)
            j = i - 1
            do while (j >= 1 .and. vertices(j) > key)
                vertices(j + 1) = vertices(j)
                j = j - 1
            end do
            vertices(j + 1) = key
        end do
    end subroutine sort_face_vertices

    integer function face_hash(vertices, hash_size) result(bucket_idx)
        integer(c_int), intent(in) :: vertices(4)
        integer, intent(in) :: hash_size
        integer(int64) :: hash
        integer :: i

        hash = 1469598103934665603_int64
        do i = 1, 4
            hash = ieor(hash, int(vertices(i), int64))
            hash = hash * 1099511628211_int64
        end do
        bucket_idx = int(mod(abs(hash), int(hash_size, int64))) + 1
    end function face_hash

end program fortran_cgal_demo
