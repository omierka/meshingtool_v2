module cleanup_mod
  use var_mod, only: element_patch_type, element_patch_group, clean_element_patches, &
                     clear_patch, clear_patch_group, rk, target_mesh, refined_clean_mesh, clear_mesh, &
                     mesh_type, refinement_depth
  use inout_mod, only: write_patch_group_vtu
  implicit none

  real(rk), parameter :: tolerance_fallback = 1.0e-10_rk
  real(rk) :: coord_tolerance = tolerance_fallback
  integer, parameter :: hex_edge_pairs(2, 12) = reshape([ &
       1, 2, 2, 3, 3, 4, 4, 1, &
       1, 5, 2, 6, 3, 7, 4, 8, &
       5, 6, 6, 7, 7, 8, 8, 5 ], [2, 12])
  type :: vertex_map_type
     integer, allocatable :: ids(:)
  end type vertex_map_type

contains

  subroutine clear_intra_patches(patches)
    type(element_patch_group), intent(in) :: patches(:)
    integer :: i

    call reset_clean_storage()
    if (size(patches) <= 0) return

    allocate(clean_element_patches(size(patches)))
    do i = 1, size(clean_element_patches)
       call clear_patch_group(clean_element_patches(i))
       call build_clean_group(patches(i), clean_element_patches(i), i)
    end do
  end subroutine clear_intra_patches
  subroutine clear_inter_patches(span_mesh)
    type(mesh_type), intent(in), optional, target :: span_mesh
    type(vertex_map_type), allocatable :: vertex_maps(:)
    real(rk), allocatable :: coords_buffer(:,:)
    integer, allocatable :: conn_buffer(:,:)
    integer, allocatable :: monitor_buffer(:)
    integer, allocatable :: candidate_ids(:)
    integer, allocatable :: neighbor_elements(:)
    integer :: total_vertices, total_cells
    integer :: elem_idx, n_elements
    integer :: global_vertex_count, cell_offset
    integer :: vertex_idx, neighbor_count, neighbor_elem
    integer :: candidate_count, local_idx, corner
    integer :: max_neighbors, i
    type(mesh_type), pointer :: adjacency_mesh
    real(rk) :: point(3)

    call clear_mesh(refined_clean_mesh)

    if (.not.allocated(clean_element_patches)) return

    if (present(span_mesh)) then
       adjacency_mesh => span_mesh
    else
       adjacency_mesh => target_mesh
    end if

    n_elements = size(clean_element_patches)
    if (n_elements <= 0) return

    total_vertices = 0
    total_cells = 0
    do elem_idx = 1, n_elements
       if (.not.allocated(clean_element_patches(elem_idx)%patchlist)) cycle
       if (clean_element_patches(elem_idx)%count <= 0) cycle
       associate(patch => clean_element_patches(elem_idx)%patchlist(1))
          if (patch%n_vert > 0) total_vertices = total_vertices + patch%n_vert
          if (patch%n_elem > 0) total_cells = total_cells + patch%n_elem
       end associate
    end do

    if (total_vertices <= 0 .or. total_cells <= 0) return

    allocate(coords_buffer(3, total_vertices))
    coords_buffer = 0.0_rk
    allocate(conn_buffer(8, total_cells))
    conn_buffer = 0
    allocate(monitor_buffer(total_cells))
    monitor_buffer = 0
    allocate(candidate_ids(total_vertices))
    candidate_ids = 0
    allocate(neighbor_elements(max(1, n_elements)))
    neighbor_elements = 0
    allocate(vertex_maps(n_elements))

    do elem_idx = 1, n_elements
       if (allocated(vertex_maps(elem_idx)%ids)) deallocate(vertex_maps(elem_idx)%ids)
       if (.not.allocated(clean_element_patches(elem_idx)%patchlist)) cycle
       if (clean_element_patches(elem_idx)%count <= 0) cycle
       associate(patch => clean_element_patches(elem_idx)%patchlist(1))
          if (patch%n_vert <= 0) cycle
          allocate(vertex_maps(elem_idx)%ids(patch%n_vert))
          vertex_maps(elem_idx)%ids = 0
       end associate
    end do

    global_vertex_count = 0
    cell_offset = 0

    do elem_idx = 1, n_elements
       if (.not.allocated(clean_element_patches(elem_idx)%patchlist)) cycle
       if (clean_element_patches(elem_idx)%count <= 0) cycle
       associate(patch => clean_element_patches(elem_idx)%patchlist(1))
          if (patch%n_vert <= 0 .or. patch%n_elem <= 0) cycle
          if (.not.allocated(patch%global_coor)) cycle
          if (.not.allocated(patch%kvert)) cycle

          neighbor_elements = 0
          neighbor_count = 1
          neighbor_elements(1) = elem_idx
          if (allocated(adjacency_mesh%kelementspan)) then
             if (elem_idx <= size(adjacency_mesh%kelementspan)) then
                if (allocated(adjacency_mesh%kelementspan(elem_idx)%list)) then
                   max_neighbors = min(size(adjacency_mesh%kelementspan(elem_idx)%list), size(neighbor_elements) - 1)
                   do i = 1, max_neighbors
                      neighbor_elements(1 + i) = adjacency_mesh%kelementspan(elem_idx)%list(i)
                   end do
                   neighbor_count = 1 + max_neighbors
                end if
             end if
          else
             neighbor_count = min(elem_idx, size(neighbor_elements))
             do i = 1, neighbor_count
                neighbor_elements(i) = elem_idx - i + 1
             end do
          end if

          do vertex_idx = 1, patch%n_vert
             candidate_count = 0
             do i = 1, neighbor_count
                neighbor_elem = neighbor_elements(i)
                if (neighbor_elem < 1 .or. neighbor_elem > n_elements) cycle
                if (.not.allocated(vertex_maps(neighbor_elem)%ids)) cycle
                candidate_count = gather_candidate_ids(vertex_maps(neighbor_elem)%ids, candidate_ids, candidate_count)
             end do
             point = patch%global_coor(:, vertex_idx)
             local_idx = find_matching_global_vertex(point, coords_buffer, candidate_ids, candidate_count, global_vertex_count)
             if (local_idx == 0) then
                global_vertex_count = global_vertex_count + 1
                coords_buffer(:, global_vertex_count) = point
                local_idx = global_vertex_count
             end if
             vertex_maps(elem_idx)%ids(vertex_idx) = local_idx
          end do

          do vertex_idx = 1, patch%n_elem
             cell_offset = cell_offset + 1
             if (allocated(patch%monitor)) then
                if (vertex_idx <= size(patch%monitor)) then
                   monitor_buffer(cell_offset) = patch%monitor(vertex_idx)
                else
                   monitor_buffer(cell_offset) = 0
                end if
             else
                monitor_buffer(cell_offset) = 0
             end if
             do corner = 1, 8
                if (corner > size(patch%kvert, 1)) then
                   conn_buffer(corner, cell_offset) = 0
                else
                   local_idx = patch%kvert(corner, vertex_idx)
                   if (local_idx >= 1 .and. local_idx <= size(vertex_maps(elem_idx)%ids)) then
                      conn_buffer(corner, cell_offset) = vertex_maps(elem_idx)%ids(local_idx)
                   else
                      conn_buffer(corner, cell_offset) = 0
                   end if
                end if
             end do
          end do
       end associate
    end do

    if (global_vertex_count > 0 .and. cell_offset > 0) then
       refined_clean_mesh%nvt = global_vertex_count
       refined_clean_mesh%nel = cell_offset
       refined_clean_mesh%nbct = 0
       refined_clean_mesh%nve = 0
       refined_clean_mesh%nee = 0
       refined_clean_mesh%nae = 0
       refined_clean_mesh%template_id = 0
       allocate(refined_clean_mesh%coor(3, global_vertex_count))
       refined_clean_mesh%coor = coords_buffer(:, 1:global_vertex_count)
       allocate(refined_clean_mesh%kvert(8, cell_offset))
       refined_clean_mesh%kvert = conn_buffer(:, 1:cell_offset)
       if (allocated(refined_clean_mesh%monitor)) deallocate(refined_clean_mesh%monitor)
       allocate(refined_clean_mesh%monitor(cell_offset))
       refined_clean_mesh%monitor = monitor_buffer(1:cell_offset)
       if (allocated(refined_clean_mesh%knpr)) deallocate(refined_clean_mesh%knpr)
       allocate(refined_clean_mesh%knpr(global_vertex_count))
       refined_clean_mesh%knpr = 0
       write(*, '(A,I0,A,I0,A,I0)') 'Inter-patch cleanup: vertices before=', total_vertices, &
            ', after=', refined_clean_mesh%nvt, ', elements=', refined_clean_mesh%nel
    else
       call clear_mesh(refined_clean_mesh)
    end if

    do elem_idx = 1, n_elements
       if (allocated(vertex_maps(elem_idx)%ids)) deallocate(vertex_maps(elem_idx)%ids)
    end do
    deallocate(vertex_maps)
    deallocate(coords_buffer, conn_buffer, monitor_buffer, candidate_ids, neighbor_elements)
  end subroutine clear_inter_patches

  subroutine reset_clean_storage()
    integer :: i

    if (.not.allocated(clean_element_patches)) return
    do i = 1, size(clean_element_patches)
       call clear_patch_group(clean_element_patches(i))
    end do
    deallocate(clean_element_patches)
  end subroutine reset_clean_storage

  subroutine build_clean_group(src_group, dst_group, root_id)
    type(element_patch_group), intent(in) :: src_group
    type(element_patch_group), intent(inout) :: dst_group
    integer, intent(in) :: root_id
    type(element_patch_type) :: merged_patch
    real(rk), allocatable :: unique_coords(:,:)
    integer, allocatable   :: unique_knpr(:)
    integer, allocatable   :: node_map(:)
    integer :: total_vertices, total_cells, unique_count
    integer :: before_count, after_count
    integer :: patch_idx, node, elem, corner, local_idx
    integer :: knpr_val, cell_offset, idx, max_depth
    real(rk) :: point(3)

    call clear_patch_group(dst_group)

    before_count = 0
    after_count = 0
    total_vertices = 0
    total_cells = 0
    max_depth = 0

    if (allocated(src_group%patchlist)) then
       do patch_idx = 1, src_group%count
          associate(src_patch => src_group%patchlist(patch_idx))
             if (src_patch%n_vert > 0) then
                total_vertices = total_vertices + src_patch%n_vert
                before_count = before_count + src_patch%n_vert
             end if
             if (src_patch%n_elem > 0) then
                total_cells = total_cells + src_patch%n_elem
             end if
             if (src_patch%depth > max_depth) max_depth = src_patch%depth
          end associate
       end do
    end if

    if (total_vertices <= 0 .or. total_cells <= 0) then
       call log_cleanup_stats(root_id, before_count, after_count)
       return
    end if

    allocate(unique_coords(3, total_vertices))
    allocate(unique_knpr(total_vertices))
    unique_coords = 0.0_rk
    unique_knpr = 0
    unique_count = 0
    cell_offset = 0

    call clear_patch(merged_patch)
    merged_patch%root_element_id = root_id
    merged_patch%template_id = 0
    merged_patch%depth = max_depth
    merged_patch%is_final = .true.
    merged_patch%vertex_pattern = .false.
    merged_patch%element_vertex_ids = 0
    merged_patch%n_elem = total_cells
    allocate(merged_patch%kvert(8, total_cells))
    merged_patch%kvert = 0
    allocate(merged_patch%monitor(total_cells))
    merged_patch%monitor = 0

    do patch_idx = 1, src_group%count
       if (.not.allocated(src_group%patchlist)) exit
       associate(src_patch => src_group%patchlist(patch_idx))
          if (src_patch%n_vert <= 0) cycle
          if (.not.allocated(src_patch%global_coor)) cycle
          allocate(node_map(src_patch%n_vert))
          do node = 1, src_patch%n_vert
             point = src_patch%global_coor(:, node)
             idx = locate_existing_vertex(point, unique_coords, unique_count)
             knpr_val = 0
             if (allocated(src_patch%knpr)) then
                if (node <= size(src_patch%knpr)) knpr_val = src_patch%knpr(node)
             end if
             if (idx == 0) then
                unique_count = unique_count + 1
                unique_coords(:, unique_count) = point
                unique_knpr(unique_count) = knpr_val
                idx = unique_count
             else
                unique_knpr(idx) = max(unique_knpr(idx), knpr_val)
             end if
             node_map(node) = idx
          end do

          if (src_patch%n_elem > 0 .and. allocated(src_patch%kvert)) then
             do elem = 1, src_patch%n_elem
                cell_offset = cell_offset + 1
                if (allocated(src_patch%monitor)) then
                   if (elem <= size(src_patch%monitor)) then
                      merged_patch%monitor(cell_offset) = src_patch%monitor(elem)
                   else
                      merged_patch%monitor(cell_offset) = 0
                   end if
                else
                   merged_patch%monitor(cell_offset) = 0
                end if
                do corner = 1, min(8, size(src_patch%kvert, 1))
                   local_idx = src_patch%kvert(corner, elem)
                   if (local_idx >= 1 .and. local_idx <= size(node_map)) then
                      merged_patch%kvert(corner, cell_offset) = node_map(local_idx)
                   else
                      merged_patch%kvert(corner, cell_offset) = 0
                   end if
                end do
             end do
          end if

          deallocate(node_map)
       end associate
    end do

    after_count = unique_count
    merged_patch%n_vert = unique_count
    if (unique_count > 0) then
       allocate(merged_patch%global_coor(3, unique_count))
       merged_patch%global_coor(:, :) = unique_coords(:, 1:unique_count)
       allocate(merged_patch%knpr(unique_count))
       merged_patch%knpr = unique_knpr(1:unique_count)
    end if

    if (cell_offset < total_cells) then
       merged_patch%n_elem = cell_offset
    end if

    call log_cleanup_stats(root_id, before_count, after_count)

    allocate(dst_group%patchlist(1))
    dst_group%patchlist(1) = merged_patch
    dst_group%count = 1

    call clear_patch(merged_patch)
    deallocate(unique_coords, unique_knpr)
  end subroutine build_clean_group

  integer function locate_existing_vertex(point, coords, count) result(idx)
    real(rk), intent(in) :: point(3)
    real(rk), intent(in) :: coords(:, :)
    integer, intent(in) :: count
    integer :: candidate
    real(rk) :: diff(3)

    idx = 0
    if (count <= 0) return
    do candidate = 1, count
       diff = coords(:, candidate) - point
       if (maxval(abs(diff)) <= coord_tolerance) then
          idx = candidate
          return
       end if
    end do
  end function locate_existing_vertex

  integer function gather_candidate_ids(source_ids, destination, start_count) result(new_count)
    integer, intent(in) :: source_ids(:)
    integer, intent(inout) :: destination(:)
    integer, intent(in) :: start_count
    integer :: i

    new_count = start_count
    do i = 1, size(source_ids)
       if (source_ids(i) <= 0) cycle
       if (new_count >= size(destination)) exit
       new_count = new_count + 1
       destination(new_count) = source_ids(i)
    end do
  end function gather_candidate_ids

  integer function find_matching_global_vertex(point, coords, candidate_ids, candidate_count, max_index) result(idx)
    real(rk), intent(in) :: point(3)
    real(rk), intent(in) :: coords(:, :)
    integer, intent(in) :: candidate_ids(:)
    integer, intent(in) :: candidate_count
    integer, intent(in) :: max_index
    integer :: i, gid
    real(rk) :: diff(3)

    idx = 0
    if (candidate_count <= 0) return
    do i = 1, candidate_count
       gid = candidate_ids(i)
       if (gid < 1 .or. gid > max_index) cycle
       diff = coords(:, gid) - point
       if (maxval(abs(diff)) <= coord_tolerance) then
          idx = gid
          return
       end if
    end do
  end function find_matching_global_vertex

  subroutine update_coord_tolerance_from_patches(patches)
    type(element_patch_group), intent(in) :: patches(:)
    real(rk) :: min_edge, edge_len, vec(3)
    integer :: grp_idx, patch_idx, elem, edge_id
    integer :: node_a, node_b

    min_edge = huge(0.0_rk)
    if (size(patches) <= 0) then
       coord_tolerance = tolerance_fallback
       return
    end if

    do grp_idx = 1, size(patches)
       if (.not.allocated(patches(grp_idx)%patchlist)) cycle
       do patch_idx = 1, patches(grp_idx)%count
          associate(patch => patches(grp_idx)%patchlist(patch_idx))
             if (patch%n_elem <= 0) cycle
             if (.not.allocated(patch%kvert)) cycle
             if (.not.allocated(patch%global_coor)) cycle
             do elem = 1, patch%n_elem
                do edge_id = 1, size(hex_edge_pairs, 2)
                   node_a = patch%kvert(hex_edge_pairs(1, edge_id), elem)
                   node_b = patch%kvert(hex_edge_pairs(2, edge_id), elem)
                   if (node_a <= 0 .or. node_a > size(patch%global_coor, 2)) cycle
                   if (node_b <= 0 .or. node_b > size(patch%global_coor, 2)) cycle
                   vec = patch%global_coor(:, node_a) - patch%global_coor(:, node_b)
                   edge_len = sqrt(sum(vec * vec))
                   if (edge_len > 0.0_rk) then
                      if (edge_len < min_edge) min_edge = edge_len
                   end if
                end do
             end do
          end associate
       end do
    end do

    if (min_edge < huge(0.0_rk) / 2.0_rk) then
      coord_tolerance = 0.5_rk * min_edge
    else
      coord_tolerance = tolerance_fallback
    end if
    write(*, '(A,ES12.5)') 'Cleanup tolerance set to ', coord_tolerance
  end subroutine update_coord_tolerance_from_patches

  subroutine write_refined_clean_vtu(mesh, filename)
    type(mesh_type), intent(in) :: mesh
    character(len=*), intent(in) :: filename
    type(element_patch_group) :: temp_group(1)
    type(element_patch_type) :: temp_patch

    if (.not.allocated(mesh%coor)) then
       write(*, '(A)') 'Refined mesh coordinates missing; VTU not written.'
       return
    end if
    if (.not.allocated(mesh%kvert)) then
       write(*, '(A)') 'Refined mesh connectivity missing; VTU not written.'
       return
    end if
    if (mesh%nvt <= 0 .or. mesh%nel <= 0) then
       write(*, '(A)') 'Refined mesh empty; VTU not written.'
       return
    end if

    call clear_patch(temp_patch)
    temp_patch%template_id = 0
    temp_patch%is_final = .true.
    temp_patch%n_vert = mesh%nvt
    temp_patch%n_elem = mesh%nel
    allocate(temp_patch%global_coor(3, temp_patch%n_vert))
    temp_patch%global_coor = mesh%coor
    allocate(temp_patch%kvert(8, temp_patch%n_elem))
    temp_patch%kvert = mesh%kvert
    allocate(temp_patch%knpr(temp_patch%n_vert))
    temp_patch%knpr = 0
    if (allocated(mesh%monitor)) then
       allocate(temp_patch%monitor(temp_patch%n_elem))
       temp_patch%monitor = mesh%monitor
    else
       allocate(temp_patch%monitor(temp_patch%n_elem))
       temp_patch%monitor = 0
    end if

    allocate(temp_group(1)%patchlist(1))
    temp_group(1)%patchlist(1) = temp_patch
    temp_group(1)%count = 1
    call write_patch_group_vtu(temp_group, filename)
    call clear_patch(temp_group(1)%patchlist(1))
    deallocate(temp_group(1)%patchlist)
  end subroutine write_refined_clean_vtu

  subroutine log_cleanup_stats(root_id, before_count, after_count)
    integer, intent(in) :: root_id
    integer, intent(in) :: before_count
    integer, intent(in) :: after_count

    return
    if (before_count > after_count) then
    write(*, '(A,I0,A,I0,A,I0)') 'Cleanup root ', root_id, ': vertices before=', &
         before_count, ', after=', after_count
    end if
  end subroutine log_cleanup_stats

end module cleanup_mod
