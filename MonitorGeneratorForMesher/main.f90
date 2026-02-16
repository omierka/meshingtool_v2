program demo
  use mpi
  use tri_tet_intersection
  use hex_io
  use cgal_bindings, only: cgal_initialize_triangle_tree, cgal_finalize_triangle_tree, cgal_have_triangle_tree
  implicit none
  integer, parameter :: max_vertices = 32
  integer, parameter :: num_tets = num_hex_tets
  real(dp) :: Hex(3,8)
  real(dp) :: subdiv_pts(3, num_hex_sub_vertices)
  real(dp) :: eps
  integer :: t, i, j, elem
  real(dp), allocatable :: dcorvg_hex(:, :)
  real(dp), allocatable :: dcorvg_tri(:, :)
  integer, allocatable :: kvert_hex(:, :)
  integer, allocatable :: kvert_tri(:, :)
  integer :: nvt_hex, nel_hex
  integer :: nvt_tri, ntri
  real(dp), allocatable :: MIS_diameter(:)
  integer :: total_tets, total_sub_vertices
  integer, allocatable :: tet_conn(:, :)
  real(dp), allocatable :: tet_points(:, :)
  integer :: idx
  integer :: total_int_cells, total_int_points, total_conn
  real(dp), allocatable :: int_points(:, :)
  integer, allocatable :: int_conn(:), int_offsets(:), int_types(:)
  integer, allocatable :: int_hex_ids(:)
  integer :: current_point, current_conn, cellIdx, nconn, start_idx
  integer :: tet_nodes(4, num_tets)
  integer :: local_vertex_indices(num_hex_sub_vertices)
  integer :: itype_tmp(num_tets), nX_tmp(num_tets)
  real(dp) :: X_tmp(3, max_vertices, num_tets)
  real(dp), allocatable :: sphere_centers(:, :)
  real(dp), allocatable :: sphere_radii(:)
  type :: triangle_index_list
    integer, allocatable :: values(:)
  end type triangle_index_list
  type(triangle_index_list), allocatable :: triangle_indices(:)
  integer, allocatable :: temp_tri_indices(:)
  integer, allocatable :: local_elem_indices(:)
  real(dp), allocatable :: hex_mis_local(:)
  real(dp), allocatable :: hex_surface_local(:)
  integer :: elem_min, elem_max, n_active, range_idx
  integer :: local_elem_min, local_elem_max, local_n_active
  integer :: tri_pos, tri_id, num_triangles_elem
  real(dp) :: triA(3), triB(3), triC(3)
  real(dp) :: mis_val
  real(dp) :: area_val
  real(dp) :: elem_mis, elem_surface,dMinGap
  logical :: element_reported
  logical :: triangle_reported
  logical :: stage_finished
  character(len=512) :: hex_file, tri_file, arg, cMinGap, output_folder
  integer :: nargs, argi
  logical :: verbose
  integer :: point_offset
  integer :: mpi_rank, mpi_size, ierr
  integer :: base_count, extra, start_offset
  integer :: local_pre_done, local_count_done, local_out_done
  character(len=512) :: file_base, tet_file_path, int_file_path
  character(len=512) :: mesh_piece_file, pvtu_path, prefix_path, area_file
  integer :: next_pct_pre, next_pct_count, next_pct_out
  logical :: progress_line_open
  real(dp) :: t_pre_start, t_pre_elapsed, t_pre_max
  real(dp) :: t_count_start, t_count_elapsed, t_count_max
  real(dp) :: t_out_start, t_out_elapsed, t_out_max
  type tMonitor
   real(dp), allocatable :: MIS(:),AREA(:),aux(:)
  end type
  type(tMonitor) :: Monitor
  integer :: ierr_comm

  eps = 1.0d-12

  verbose = .false.
  hex_file = "INPUT/hex.tri"
  tri_file = "INPUT/tri.off"
  output_folder = "OUTPUT"
  dMinGap = -1.0_dp
  nargs = command_argument_count()
  argi = 1
  do while (argi <= nargs)
    call get_command_argument(argi, arg)
    select case (trim(arg))
    case ('-h','--hex')
      if (argi + 1 > nargs) then
        print *, 'Missing value for ', trim(arg)
        stop 1
      end if
      call get_command_argument(argi+1, hex_file)
      argi = argi + 1
    case ('-t','--tri')
      if (argi + 1 > nargs) then
        print *, 'Missing value for ', trim(arg)
        stop 1
      end if
      call get_command_argument(argi+1, tri_file)
      argi = argi + 1
    case ('-o','--output-folder')
      if (argi + 1 > nargs) then
        print *, 'Missing value for ', trim(arg)
        stop 1
      end if
      call get_command_argument(argi+1, output_folder)
      argi = argi + 1
    case ('-m','--mingap')
      if (argi + 1 > nargs) then
        print *, 'Missing value for ', trim(arg)
        stop 1
      end if
      call get_command_argument(argi+1, cMinGap)
      Read(cMinGap,*) dMinGap
      argi = argi + 1
    case ('-v','--verbose')
      verbose = .true.
    case default
      print *, 'Unknown argument: ', trim(arg)
      stop 1
    end select
    argi = argi + 1
  end do

  progress_line_open = .false.

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)
  if (mpi_rank == 0) call ensure_directory(trim(output_folder))
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  call read_hex_mesh(trim(hex_file), dcorvg_hex, kvert_hex, nvt_hex, nel_hex)
  if (mpi_rank == 0) then
    write(*,'(A,1X,A,1X,"with",1X,I0,1X,"vertices and",1X,I0,1X,"hexahedra")') &
         'Loaded hex mesh:', trim(hex_file), nvt_hex, nel_hex
  end if

  call read_tri_mesh(trim(tri_file), dcorvg_tri, kvert_tri, nvt_tri, ntri, MIS_diameter)
  do i = 1, size(MIS_diameter)
    if (MIS_diameter(i) < 0.0_dp) MIS_diameter(i) = 10.0_dp*dMinGap
  end do
  if (mpi_rank == 0) then
    write(*,'(A,1X,A,1X,"with",1X,I0,1X,"vertices and",1X,I0,1X,"triangles")') &
         'Loaded triangle mesh:', trim(tri_file), nvt_tri, ntri
  end if

  call cgal_initialize_triangle_tree(nvt_tri, ntri, dcorvg_tri, kvert_tri)
  if (mpi_rank == 0) then
    if (cgal_have_triangle_tree()) then
      write(*,*) 'CGAL triangle accelerator initialized (AABB tree ready)'
    else
      write(*,*) 'CGAL accelerator unavailable'
    end if
  end if
  if (.not. cgal_have_triangle_tree()) then
    if (mpi_rank == 0) write(*,*) 'CGAL acceleration required for this build; aborting.'
    call cgal_finalize_triangle_tree()
    call MPI_Finalize(ierr)
    stop 1
  end if

!   elem_min = 37
!   elem_max = 37
  elem_min = 1
  elem_max = nel_hex

  if (dMinGap < 0.0_dp) then
    if (mpi_rank == 0) print *, "the --mingap/-m parameter for mingap is not properly set"
    call cgal_finalize_triangle_tree()
    call MPI_Finalize(ierr)
    stop 1
  end if

  if (nel_hex <= 0) then
    if (mpi_rank == 0) print *, "Mesh contains no hexahedra"
    call cgal_finalize_triangle_tree()
    call MPI_Finalize(ierr)
    stop 1
  end if

  allocate(Monitor%MIS(nel_hex),Monitor%AREA(nel_hex),Monitor%aux(nel_hex))
  Monitor%MIS = 0.0_dp
  Monitor%AREA = 0.0_dp

  if (elem_min < 1) elem_min = 1
  if (elem_max > nel_hex) elem_max = nel_hex
  if (elem_min > elem_max) then
    if (mpi_rank == 0) print *, "Invalid element range"
    call cgal_finalize_triangle_tree()
    call MPI_Finalize(ierr)
    stop 1
  end if
  n_active = elem_max - elem_min + 1

  base_count = n_active / mpi_size
  extra = mod(n_active, mpi_size)
  start_offset = mpi_rank * base_count + min(mpi_rank, extra)
  local_n_active = base_count
  if (mpi_rank < extra) local_n_active = local_n_active + 1
  local_elem_min = elem_min + start_offset
  if (local_n_active > 0) then
    local_elem_max = local_elem_min + local_n_active - 1
  else
    local_elem_max = local_elem_min - 1
  end if

  next_pct_pre = 1
  next_pct_count = 1
  next_pct_out = 1

  total_tets = num_tets * local_n_active
  total_sub_vertices = num_hex_sub_vertices * local_n_active
  allocate(tet_conn(4, total_tets))
  if (total_sub_vertices > 0) then
    allocate(tet_points(3, total_sub_vertices))
  else
    allocate(tet_points(3,0))
  end if
  if (local_n_active > 0) then
    allocate(sphere_centers(3, local_n_active))
    allocate(sphere_radii(local_n_active))
    allocate(triangle_indices(local_n_active))
    allocate(hex_mis_local(local_n_active))
    allocate(hex_surface_local(local_n_active))
    hex_mis_local = 0.0_dp
    hex_surface_local = 0.0_dp
  else
    allocate(sphere_centers(3,0))
    allocate(sphere_radii(0))
    allocate(triangle_indices(0))
    allocate(hex_mis_local(0))
    allocate(hex_surface_local(0))
  end if

  local_pre_done = 0
  stage_finished = (n_active == 0)
  if (.not. stage_finished) then
    if (.not. verbose .and. mpi_rank == 0) write(*,*) 'Preprocessing hexahedra...'
    elem = local_elem_min
  end if
  t_pre_start = MPI_Wtime()
  do while (.not. stage_finished)
    if (elem <= local_elem_max) then
      range_idx = elem - local_elem_min + 1
      do j = 1, 8
        Hex(:, j) = dcorvg_hex(:, kvert_hex(j, elem))
      end do
      call compute_hex_subdivision_points(Hex, subdiv_pts)
      point_offset = (range_idx-1)*num_hex_sub_vertices
      if (total_sub_vertices > 0) then
        tet_points(:, point_offset+1:point_offset+num_hex_sub_vertices) = subdiv_pts
      end if
      do j = 1, num_hex_sub_vertices
        local_vertex_indices(j) = point_offset + j
      end do
      call compute_hex_bounding_sphere(Hex, sphere_centers(:,range_idx), sphere_radii(range_idx))
      call gather_triangles_in_sphere(sphere_centers(:,range_idx), sphere_radii(range_idx), Hex, dcorvg_tri, &
           kvert_tri, ntri, temp_tri_indices)
      if (allocated(triangle_indices(range_idx)%values)) deallocate(triangle_indices(range_idx)%values)
      triangle_indices(range_idx)%values = temp_tri_indices
      if (verbose) write(*,*) 'Element', elem, 'candidate triangles:', temp_tri_indices
      if (allocated(temp_tri_indices)) deallocate(temp_tri_indices)
      call hexa_tet_node_indices(local_vertex_indices, tet_nodes)
      do t = 1, num_tets
        idx = (range_idx-1)*num_tets + t
        tet_conn(:, idx) = tet_nodes(:, t)
      end do
      local_pre_done = local_pre_done + 1
      elem = elem + 1
    end if
    call sync_progress(local_pre_done, n_active, next_pct_pre, '  Preprocessing progress: ', mpi_rank, stage_finished)
  end do
  t_pre_elapsed = MPI_Wtime() - t_pre_start
  call MPI_Reduce(t_pre_elapsed, t_pre_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  if (mpi_rank == 0) write(*,'(A,F10.3," s")') 'Preprocessing loop time (max rank):', t_pre_max

  total_int_cells = 0
  total_int_points = 0
  total_conn = 0
  local_count_done = 0
  stage_finished = (n_active == 0)
  if (.not. stage_finished) then
    if (.not. verbose .and. mpi_rank == 0) write(*,*) 'Counting potential intersections...'
    elem = local_elem_min
  end if
  t_count_start = MPI_Wtime()
  do while (.not. stage_finished)
    if (elem <= local_elem_max) then
      range_idx = elem - local_elem_min + 1
      do j = 1, 8
        Hex(:, j) = dcorvg_hex(:, kvert_hex(j, elem))
      end do
      if (allocated(triangle_indices(range_idx)%values)) then
        num_triangles_elem = size(triangle_indices(range_idx)%values)
        if (num_triangles_elem > 0) then
          do tri_pos = 1, num_triangles_elem
            tri_id = triangle_indices(range_idx)%values(tri_pos)
            triA = dcorvg_tri(:, kvert_tri(1, tri_id))
            triB = dcorvg_tri(:, kvert_tri(2, tri_id))
            triC = dcorvg_tri(:, kvert_tri(3, tri_id))
            call tri_hex_intersections(triA, triB, triC, Hex, eps, itype_tmp, nX_tmp, X_tmp)
            do t = 1, num_tets
              if (nX_tmp(t) <= 0 .or. itype_tmp(t) <= 0) cycle
              total_int_cells = total_int_cells + 1
              total_int_points = total_int_points + nX_tmp(t)
              select case (itype_tmp(t))
              case (1)
                total_conn = total_conn + 1
              case (2)
                total_conn = total_conn + 2
              case default
                total_conn = total_conn + nX_tmp(t)
              end select
            end do
          end do
        end if
      end if
      local_count_done = local_count_done + 1
      elem = elem + 1
    end if
    call sync_progress(local_count_done, n_active, next_pct_count, '  Counting progress: ', mpi_rank, stage_finished)
  end do
  t_count_elapsed = MPI_Wtime() - t_count_start
  call MPI_Reduce(t_count_elapsed, t_count_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  if (mpi_rank == 0) write(*,'(A,F10.3," s")') 'Counting loop time (max rank):   ', t_count_max

  if (total_int_points > 0) then
    allocate(int_points(3, total_int_points))
  else
    allocate(int_points(3,0))
  end if
  if (total_conn > 0) then
    allocate(int_conn(total_conn))
  else
    allocate(int_conn(0))
  end if
  if (total_int_cells > 0) then
    allocate(int_offsets(total_int_cells), int_types(total_int_cells), int_hex_ids(total_int_cells))
  else
    allocate(int_offsets(0), int_types(0), int_hex_ids(0))
  end if

  current_point = 0
  current_conn = 1
  cellIdx = 1
  stage_finished = (n_active == 0)
  if (.not. stage_finished) then
    if (.not. verbose .and. mpi_rank == 0) write(*,*) 'Computing intersections...'
    elem = local_elem_min
  end if
  local_out_done = 0
  t_out_start = MPI_Wtime()
  do while (.not. stage_finished)
    if (elem <= local_elem_max) then
      elem_surface = 0.0_dp
      elem_mis = 1d8
      range_idx = elem - local_elem_min + 1
      do j = 1, 8
        Hex(:, j) = dcorvg_hex(:, kvert_hex(j, elem))
      end do
      if (allocated(triangle_indices(range_idx)%values)) then
        num_triangles_elem = size(triangle_indices(range_idx)%values)
        if (num_triangles_elem > 0) then
          element_reported = .false.
          do tri_pos = 1, num_triangles_elem
            tri_id = triangle_indices(range_idx)%values(tri_pos)
            triA = dcorvg_tri(:, kvert_tri(1, tri_id))
            triB = dcorvg_tri(:, kvert_tri(2, tri_id))
            triC = dcorvg_tri(:, kvert_tri(3, tri_id))
            call tri_hex_intersections(triA, triB, triC, Hex, eps, itype_tmp, nX_tmp, X_tmp)
            triangle_reported = .false.
            do t = 1, num_tets
              if (nX_tmp(t) <= 0 .or. itype_tmp(t) <= 0) cycle
              if (.not. element_reported .and. verbose) then
                print *, "Element", elem
                print *, "  Sphere center=", sphere_centers(:,range_idx), " radius=", sphere_radii(range_idx)
                print *, "  Triangles inside sphere:", num_triangles_elem
                element_reported = .true.
              end if
              if (.not. triangle_reported .and. verbose) then
                print *, "  Triangle", tri_id
                triangle_reported = .true.
              end if
              if (verbose) print *, "    Tetra", t, "itype=", itype_tmp(t), " nX=", nX_tmp(t)
              do i = 1, nX_tmp(t)
                if (verbose) print *, "      ", i, X_tmp(:,i,t)
                current_point = current_point + 1
                int_points(:, current_point) = X_tmp(:, i, t)
              end do
              start_idx = current_point - nX_tmp(t) + 1
              select case (itype_tmp(t))
              case (1)
                nconn = min(1, nX_tmp(t))
                int_types(cellIdx) = 1
              case (2)
                nconn = min(2, nX_tmp(t))
                int_types(cellIdx) = 3
              case default
                nconn = nX_tmp(t)
                int_types(cellIdx) = 7
                mis_val = MIS_diameter(tri_id)
                if (mis_val < elem_mis) elem_mis = mis_val
                area_val = polygon_area(X_tmp(:,:,t), nX_tmp(t))
                elem_surface = elem_surface + area_val
              end select
              int_hex_ids(cellIdx) = elem
              do j = 1, nconn
                int_conn(current_conn) = start_idx + j - 2
                current_conn = current_conn + 1
              end do
              int_offsets(cellIdx) = current_conn - 1
              cellIdx = cellIdx + 1
            end do
          end do
        end if
      end if
      if (range_idx >= 1 .and. range_idx <= local_n_active) then
        if (elem_mis < 0.0_dp) then
          hex_mis_local(range_idx) = 0.0_dp
          Monitor%MIS(local_elem_min + range_idx - 1) = 0.0_dp
        else
          hex_mis_local(range_idx) = elem_mis
          Monitor%MIS(local_elem_min + range_idx - 1) = elem_mis
        end if
        Monitor%AREA(local_elem_min + range_idx - 1) = elem_surface
      end if
      local_out_done = local_out_done + 1
!      write(*,*) 'EL : ', elem, range_idx, elem_surface,elem_mis
      elem = elem + 1
    end if
    call sync_progress(local_out_done, n_active, next_pct_out, '  Intersection progress: ', mpi_rank, stage_finished)
  end do
  t_out_elapsed = MPI_Wtime() - t_out_start
  call MPI_Reduce(t_out_elapsed, t_out_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  if (mpi_rank == 0) write(*,'(A,F10.3," s")') 'Intersection loop time (max rank):', t_out_max

  call MPI_Allreduce(Monitor%MIS, Monitor%aux, nel_hex, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr_comm)
  Monitor%MIS = Monitor%aux
  call MPI_Allreduce(Monitor%AREA, Monitor%aux, nel_hex, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr_comm)
  Monitor%AREA = Monitor%aux
  if (mpi_rank == 0) then
    call compose_path(output_folder, 'area.txt', area_file)
    open(file=trim(area_file),unit=11)
    do i=1,nel_hex
      write(11,'(ES12.4)') 5.0_dp*dMinGap / Monitor%MIS(i)
    end do
    close(11)
  end if

  if (local_n_active > 0) then
    allocate(local_elem_indices(local_n_active))
    do i = 1, local_n_active
      local_elem_indices(i) = local_elem_min + i - 1
    end do
  else
    allocate(local_elem_indices(0))
  end if
  write(mesh_piece_file,'("hex_mesh_rank",I4.4,".vtu")') mpi_rank
  call compose_path(output_folder, trim(mesh_piece_file), mesh_piece_file)
  call write_hex_mesh_vtu(trim(mesh_piece_file), dcorvg_hex, kvert_hex, local_elem_indices, hex_mis_local, "MIS_diameter", &
       hex_surface_local, "SurfaceArea")
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (mpi_rank == 0) then
    call compose_path(output_folder, 'hex_mesh.pvtu', pvtu_path)
    call compose_path(output_folder, 'hex_mesh_rank', prefix_path)
    call write_pvtu_reference(trim(pvtu_path), mpi_size, trim(prefix_path), ".vtu", &
         cell_fields=(/'MIS_diameter','SurfaceArea '/))
  end if
  if (allocated(local_elem_indices)) deallocate(local_elem_indices)

  write(file_base,'("hex_intersection_rank",I4.4)') mpi_rank
  call compose_path(output_folder, trim(file_base), file_base)
  call write_hex_intersection_vtu(trim(file_base), tet_points, tet_conn, &
       int_points, int_conn, int_offsets, int_types, int_hex_ids, tet_file_path, int_file_path)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (mpi_rank == 0) then
    call compose_path(output_folder, 'hex_intersection_tets.pvtu', pvtu_path)
    call compose_path(output_folder, 'hex_intersection_rank', prefix_path)
    call write_pvtu_reference(trim(pvtu_path), mpi_size, trim(prefix_path), "_tets.vtu")
    call compose_path(output_folder, 'hex_intersection_intersections.pvtu', pvtu_path)
    call write_pvtu_reference(trim(pvtu_path), mpi_size, trim(prefix_path), "_intersections.vtu", "hex_id")
  end if

  call cgal_finalize_triangle_tree()
  call MPI_Finalize(ierr)
contains

  subroutine sync_progress(local_done, total_work, next_pct, label, rank, finished)
    integer, intent(in) :: local_done, total_work, rank
    integer, intent(inout) :: next_pct
    character(len=*), intent(in) :: label
    logical, intent(inout) :: finished
    integer :: global_done, ierr_loc, pct_loc

    if (finished) return
    if (total_work <= 0) then
      finished = .true.
      return
    end if

    call MPI_Allreduce(local_done, global_done, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr_loc)
    if (rank == 0) then
      pct_loc = int( (real(global_done, dp)/real(total_work, dp))*100.0_dp )
      pct_loc = max(0, min(100, pct_loc))
      do while (pct_loc >= next_pct .and. next_pct <= 100)
        call ensure_progress_line(label)
        write(*,'(A)', advance='no') '%'
        next_pct = next_pct + 1
      end do
      if (global_done >= total_work) call finalize_progress_line()
    end if
    if (global_done >= total_work) finished = .true.
  end subroutine sync_progress

  subroutine ensure_progress_line(label)
    character(len=*), intent(in) :: label
    if (.not. progress_line_open) then
      write(*,'(A)', advance='no') trim(label)//' ['
      progress_line_open = .true.
    end if
  end subroutine ensure_progress_line

  subroutine finalize_progress_line()
    if (progress_line_open) then
      write(*,'(A)', advance='no') '] 100%'
      write(*,*)
      progress_line_open = .false.
    end if
  end subroutine finalize_progress_line

  subroutine ensure_directory(dir)
    character(len=*), intent(in) :: dir
    integer :: status
    character(len=1024) :: cmd
    if (len_trim(dir) == 0) return
    if (len_trim(dir) == 1 .and. dir(1:1) == '.') return
    cmd = 'mkdir -p "'//trim(dir)//'"'
    call execute_command_line(cmd, exitstat=status)
    if (status /= 0) then
      write(*,*) 'Warning: failed to create output folder ', trim(dir)
    end if
  end subroutine ensure_directory

  subroutine compose_path(dir, name, out_path)
    character(len=*), intent(in) :: dir, name
    character(len=*), intent(out) :: out_path
    integer :: len_dir
    out_path = ''
    len_dir = len_trim(dir)
    if (len_dir == 0 .or. (len_dir == 1 .and. dir(1:1) == '.')) then
      out_path = trim(name)
    else
      if (dir(len_dir:len_dir) == '/' .or. dir(len_dir:len_dir) == '\') then
        out_path = trim(dir)//trim(name)
      else
        out_path = trim(dir)//'/'//trim(name)
      end if
    end if
  end subroutine compose_path

  pure function polygon_area(points, npts) result(area)
    real(dp), intent(in) :: points(:, :)
    integer, intent(in) :: npts
    real(dp) :: area
    real(dp) :: ref(3), v1(3), v2(3), cross_vec(3)
    integer :: k

    if (npts < 3) then
      area = 0.0_dp
      return
    end if
    ref = points(:,1)
    cross_vec = 0.0_dp
    do k = 2, npts-1
      v1 = points(:,k) - ref
      v2 = points(:,k+1) - ref
      cross_vec = cross_vec + cross_product(v1, v2)
    end do
    area = 0.5_dp * sqrt(sum(cross_vec**2))
  end function polygon_area

  pure function cross_product(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product

end program demo
