module def_mod
  use inout_mod, only: load_all_meshes, release_all_meshes, load_target_mesh
  use var_mod, only: templates, target_mesh, template_is_final, meshes, element_patches, &
                     element_patch_type, clear_patch, release_element_patches, mesh_file_count, rk, &
                     rotate_patch, determine_template, element_patch_count, release_element_span, &
                     clear_patch_group, append_patch_to_group, clear_mesh, refined_clean_mesh, mesh_type, &
                     release_mesh_levels, reproducibility, refinement_depth, Monitor_threshold, &
                     apply_cylindric_transform, cylindrical_outer_radius
  use preprocessor_config_mod, only: get_fullcyl_min_radius_percentage, get_fullcyl_max_radius_percentage
  use setupe3dfile_reader, only: MeshConfig, ProcessParameters, initialize_mesh_config, load_mesh_config, &
                                 initialize_process_parameters, load_process_parameters
  use bc_treatment, only: FaceList, BoxBoundaryClassification, HollowCylinderBoundaryClassification, &
                          InflowBoundaryGroup, recompute_knpr_from_connectivity, classify_box_boundaries, &
                          classify_hollowcylinder_boundaries
  use iso_c_binding, only: c_double, c_int
  implicit none
  integer, parameter :: max_refinement_depth = 2
  integer, parameter :: face_vertex_map(4, 6) = reshape([ &
       1, 2, 3, 4, &
       1, 2, 5, 4, &
       2, 3, 6, 5, &
       3, 4, 8, 7, &
       4, 1, 5, 8, &
       5, 6, 7, 8], [4, 6])

contains

  subroutine initialize_mesh_database()
    call load_all_meshes()
  end subroutine initialize_mesh_database

  subroutine finalize_mesh_database()
    call release_element_patches()
    call release_mesh_levels()
    call release_all_meshes()
  end subroutine finalize_mesh_database

  subroutine load_target_mesh_file(filename)
    character(len=*), intent(in) :: filename

    call load_target_mesh(filename)
  end subroutine load_target_mesh_file

  subroutine random_element_marking(percent)
    integer, intent(in) :: percent
    integer :: elem, seed_size, clock_count, i
    integer, allocatable :: seed(:)
    real(rk) :: threshold, rnd_val, monitor_rand

    if (percent <= 0) then
       if (allocated(target_mesh%monitor)) target_mesh%monitor = 0
       return
    end if
    if (target_mesh%nel <= 0) then
       write(*, '(A)') 'Target mesh has no elements; skipping random element marking.'
       return
    end if

    threshold = real(percent, rk) / 100.0_rk
    if (.not.allocated(target_mesh%monitor)) then
       allocate(target_mesh%monitor(target_mesh%nel))
    end if
    target_mesh%monitor = 0

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    if (reproducibility) then
       seed = [(123456 + 17 * i, i = 1, seed_size)]
    else
       call system_clock(count=clock_count)
       if (clock_count == 0) clock_count = 8675309
       do i = 1, seed_size
          seed(i) = clock_count + 37 * i
       end do
    end if
    call random_seed(put=seed)
    deallocate(seed)

    if (threshold >= 1.0_rk) then
       do elem = 1, target_mesh%nel
          call random_number(monitor_rand)
          target_mesh%monitor(elem) = merge(1, 2, monitor_rand < 0.2_rk)
       end do
       return
    end if

    do elem = 1, target_mesh%nel
       call random_number(rnd_val)
       if (rnd_val < threshold) then
          call random_number(monitor_rand)
          target_mesh%monitor(elem) = merge(1, 2, monitor_rand < 0.2_rk)
       end if
    end do
  end subroutine random_element_marking

  subroutine load_monitor_file(mesh, filename)
    type(mesh_type), intent(inout) :: mesh
    character(len=*), intent(in) :: filename
    integer :: unit, ios, elem
    real(rk) :: value
    logical :: exists
    character(len=128) :: label

    inquire(file=trim(filename), exist=exists)
    if (.not.exists) then
       write(*, '(A)') 'Monitor file not found: ' // trim(filename)
       return
    end if
    if (mesh%nel <= 0) then
       write(*, '(A)') 'Monitor file ignored: mesh has no elements.'
       return
    end if

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*, '(A,I0)') 'Unable to open monitor file, IOSTAT=', ios
       return
    end if

    if (.not.allocated(mesh%monitor)) allocate(mesh%monitor(mesh%nel))
    if (.not.allocated(mesh%monitor_value)) allocate(mesh%monitor_value(mesh%nel))
    mesh%monitor = 0
    mesh%monitor_value = 0.0_rk

    do elem = 1, mesh%nel
       read(unit, *, iostat=ios) value
       if (ios /= 0) then
          write(*, '(A,I0)') 'Insufficient monitor entries; stopping at element ', elem
          exit
       end if
       mesh%monitor_value(elem) = value
       if (value > Monitor_threshold) then
          mesh%monitor(elem) = refinement_depth
       else
          mesh%monitor(elem) = 0
       end if
    end do
    close(unit)

    write(label, '(A,I0,A)') 'No Of Elements with Refinement depth [0..', refinement_depth, ']'
    call report_refinement_distribution(mesh, trim(label))
  end subroutine load_monitor_file

  subroutine vertice_marking(mesh, min_marker)
    type(mesh_type), intent(inout) :: mesh
    integer, intent(in) :: min_marker
    integer :: elem, vert_idx, node_id
    integer :: marked_count

!    write(*,*) allocated(mesh%monitor),allocated(mesh%kvert),allocated(mesh%knpr)

    if (.not.allocated(mesh%monitor)) return
    if (.not.allocated(mesh%kvert)) return
    if (.not.allocated(mesh%knpr)) return

    mesh%knpr = 0
    marked_count = 0

    do elem = 1, size(mesh%kvert, 2)
       if (elem > size(mesh%monitor)) exit
       if (mesh%monitor(elem) < min_marker) cycle
       marked_count = marked_count + 1
       do vert_idx = 1, min(8, size(mesh%kvert, 1))
          node_id = mesh%kvert(vert_idx, elem)
          if (node_id > 0 .and. node_id <= size(mesh%knpr)) then
             mesh%knpr(node_id) = 1
          end if
       end do
    end do

    write(*, '(A,I0,A,I0)') 'Elements with monitor >= ', min_marker, ': ', marked_count
  end subroutine vertice_marking

  subroutine build_target_element_span(mesh)
    type(mesh_type), intent(inout) :: mesh
    integer :: nel, nvt, elem, vtx, node
    integer :: total_refs, offset, degree
    integer :: other_elem, neighbor_count, marker, i
    integer, allocatable :: vertex_counts(:), vertex_offsets(:), vertex_elems(:)
    integer, allocatable :: fill_counts(:), neighbor_marks(:), neighbor_buffer(:)

    if (.not.allocated(mesh%kvert)) return
    nel = mesh%nel
    nvt = mesh%nvt
    if (nel <= 0 .or. nvt <= 0) return

    allocate(vertex_counts(nvt))
    vertex_counts = 0
    total_refs = 0

    do elem = 1, nel
       do vtx = 1, min(8, size(mesh%kvert, 1))
          node = mesh%kvert(vtx, elem)
          if (node > 0 .and. node <= nvt) then
             vertex_counts(node) = vertex_counts(node) + 1
             total_refs = total_refs + 1
          end if
       end do
    end do

    if (allocated(mesh%kelementspan)) then
       call release_element_span(mesh%kelementspan)
       deallocate(mesh%kelementspan)
    end if
    allocate(mesh%kelementspan(nel))

    if (total_refs == 0) then
       deallocate(vertex_counts)
       return
    end if

    allocate(vertex_offsets(nvt + 1))
    vertex_offsets(1) = 1
    do node = 1, nvt
       vertex_offsets(node + 1) = vertex_offsets(node) + vertex_counts(node)
    end do
    allocate(vertex_elems(total_refs))
    allocate(fill_counts(nvt))
    fill_counts = 0

    do elem = 1, nel
       do vtx = 1, min(8, size(mesh%kvert, 1))
          node = mesh%kvert(vtx, elem)
          if (node > 0 .and. node <= nvt) then
             offset = vertex_offsets(node) + fill_counts(node)
             vertex_elems(offset) = elem
             fill_counts(node) = fill_counts(node) + 1
          end if
       end do
    end do

    allocate(neighbor_marks(nel))
    neighbor_marks = 0
    allocate(neighbor_buffer(max(1, nel)))
    marker = 0

    do elem = 1, nel
       marker = marker + 1
       neighbor_count = 0
       do vtx = 1, min(8, size(mesh%kvert, 1))
          node = mesh%kvert(vtx, elem)
          if (node <= 0 .or. node > nvt) cycle
          degree = vertex_counts(node)
          if (degree <= 0) cycle
          offset = vertex_offsets(node)
          do i = 0, degree - 1
             other_elem = vertex_elems(offset + i)
             if (other_elem == elem) cycle
             if (neighbor_marks(other_elem) == marker) cycle
             neighbor_marks(other_elem) = marker
             neighbor_count = neighbor_count + 1
             neighbor_buffer(neighbor_count) = other_elem
          end do
       end do
       if (neighbor_count > 0) then
          allocate(mesh%kelementspan(elem)%list(neighbor_count))
          mesh%kelementspan(elem)%list = neighbor_buffer(1:neighbor_count)
       else
          if (allocated(mesh%kelementspan(elem)%list)) deallocate(mesh%kelementspan(elem)%list)
       end if
    end do

    deallocate(vertex_counts, vertex_offsets, vertex_elems, fill_counts, neighbor_marks, neighbor_buffer)
  end subroutine build_target_element_span

  subroutine build_face_adjacency(mesh)
    type(mesh_type), intent(inout) :: mesh
    integer :: nel, total_faces, elem, face, node, face_idx
    integer :: group_start, group_end
    integer, allocatable :: face_keys(:, :)
    integer, allocatable :: face_owner(:), face_local(:), order(:)

    if (.not.allocated(mesh%kvert)) return
    nel = mesh%nel
    if (nel <= 0) return

    total_faces = 6 * nel
    if (allocated(mesh%kadj)) then
       deallocate(mesh%kadj)
    end if
    allocate(mesh%kadj(6, nel))
    mesh%kadj = 0

    allocate(face_keys(4, total_faces))
    allocate(face_owner(total_faces))
    allocate(face_local(total_faces))
    allocate(order(total_faces))

    face_idx = 0
    do elem = 1, nel
       do face = 1, 6
          face_idx = face_idx + 1
          do node = 1, 4
             face_keys(node, face_idx) = mesh%kvert(face_vertex_map(node, face), elem)
          end do
          call sort_face_nodes(face_keys(:, face_idx))
          face_owner(face_idx) = elem
          face_local(face_idx) = face
          order(face_idx) = face_idx
       end do
    end do

    call sort_face_records(face_keys, order)

    face_idx = 1
    do while (face_idx <= total_faces)
       group_start = face_idx
       group_end = face_idx
       do while (group_end < total_faces)
          if (.not.faces_have_same_key(face_keys, order(group_end), order(group_end + 1))) exit
          group_end = group_end + 1
       end do
       if (group_end == group_start + 1) then
          call connect_face_pair(face_owner(order(group_start)), face_local(order(group_start)), &
                                 face_owner(order(group_end)), face_local(order(group_end)))
       else if (group_end > group_start + 1) then
          write(*, '(A)') 'Warning: more than two faces share identical connectivity; adjacency may be ambiguous.'
          call connect_face_pair(face_owner(order(group_start)), face_local(order(group_start)), &
                                 face_owner(order(group_start + 1)), face_local(order(group_start + 1)))
       end if
       face_idx = group_end + 1
    end do

    deallocate(face_keys, face_owner, face_local, order)

  contains

    subroutine connect_face_pair(elem_a, face_a, elem_b, face_b)
      integer, intent(in) :: elem_a, face_a, elem_b, face_b

      if (elem_a < 1 .or. elem_a > size(mesh%kadj, 2)) return
      if (elem_b < 1 .or. elem_b > size(mesh%kadj, 2)) return
      if (face_a < 1 .or. face_a > size(mesh%kadj, 1)) return
      if (face_b < 1 .or. face_b > size(mesh%kadj, 1)) return

      mesh%kadj(face_a, elem_a) = elem_b
      mesh%kadj(face_b, elem_b) = elem_a
    end subroutine connect_face_pair

    subroutine sort_face_records(keys, order)
      integer, intent(in) :: keys(:, :)
      integer, intent(inout) :: order(:)
      integer :: i, j, key

      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (.not.face_key_less(keys, key, order(j))) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
    end subroutine sort_face_records

    logical function face_key_less(keys, ia, ib)
      integer, intent(in) :: keys(:, :)
      integer, intent(in) :: ia, ib
      integer :: idx

      face_key_less = .false.
      do idx = 1, size(keys, 1)
         if (keys(idx, ia) < keys(idx, ib)) then
            face_key_less = .true.
            return
         else if (keys(idx, ia) > keys(idx, ib)) then
            face_key_less = .false.
            return
         end if
      end do
    end function face_key_less

    logical function faces_have_same_key(keys, ia, ib)
      integer, intent(in) :: keys(:, :)
      integer, intent(in) :: ia, ib

      faces_have_same_key = all(keys(:, ia) == keys(:, ib))
    end function faces_have_same_key

  end subroutine build_face_adjacency

  subroutine sort_face_nodes(nodes)
    integer, intent(inout) :: nodes(4)
    integer :: i, j, tmp

    do i = 1, 3
       do j = i + 1, 4
          if (nodes(j) < nodes(i)) then
             tmp = nodes(i)
             nodes(i) = nodes(j)
             nodes(j) = tmp
          end if
       end do
    end do
  end subroutine sort_face_nodes

  subroutine mesh_refinement(mesh)
    type(mesh_type), intent(inout) :: mesh
    integer :: elem, i
    integer :: used_templates(mesh_file_count)
    integer :: unknown_templates
    integer :: ivert(8)
    integer :: monitor_value

    if (.not.allocated(mesh%kvert) .or. .not.allocated(mesh%knpr)) then
       write(*, '(A)') 'Target mesh not loaded; skipping refinement.'
       return
    end if
    if (.not.allocated(mesh%coor)) then
       write(*, '(A)') 'Target mesh coordinates missing; skipping refinement.'
       return
    end if
    if (.not.allocated(meshes)) then
       write(*, '(A)') 'Template database not initialized; skipping refinement.'
       return
    end if

    used_templates = 0
    unknown_templates = 0

    if (allocated(element_patches)) call release_element_patches()
    if (mesh%nel > 0) then
       allocate(element_patches(mesh%nel))
       do i = 1, mesh%nel
          call clear_patch_group(element_patches(i))
       end do
    else
       allocate(element_patches(0))
    end if
    element_patch_count = 0

    do elem = 1, mesh%nel
       ivert = mesh%kvert(:, elem)
       monitor_value = 0
       if (allocated(mesh%monitor)) then
          if (elem <= size(mesh%monitor)) monitor_value = mesh%monitor(elem)
       end if
       call refine_element(elem, mesh%coor, mesh%knpr, ivert, monitor_value, 0, &
            used_templates, unknown_templates)
    end do

    call write_template_usage(used_templates)
    write(*, '(A,I0)') ' unknown templates: ', unknown_templates
  end subroutine mesh_refinement

  recursive subroutine refine_element(root_element_id, coords_src, knpr_src, ivert, monitor_value, &
       depth, used_templates, unknown_templates)
    integer, intent(in) :: root_element_id
    real(rk), intent(in) :: coords_src(:, :)
    integer, intent(in) :: knpr_src(:)
    integer, intent(in) :: ivert(8)
    integer, intent(in) :: monitor_value
    integer, intent(in) :: depth
    integer, intent(inout) :: used_templates(:)
    integer, intent(inout) :: unknown_templates
    logical :: bvert(8), work_c(8)
    integer :: work_ids(8)
    integer :: template_id, i
    logical :: bfinal
    type(element_patch_type) :: patch

    do i = 1, 8
       if (ivert(i) > 0 .and. ivert(i) <= size(knpr_src)) then
          bvert(i) = knpr_src(ivert(i)) /= 0
       else
          bvert(i) = .false.
       end if
    end do

    work_c = bvert
    work_ids = ivert

    call rotate_patch(work_c, work_ids)
    call determine_template(work_c, template_id, bfinal)

    if (template_id <= 0) then
       unknown_templates = unknown_templates + 1
       return
    end if

    used_templates(template_id) = used_templates(template_id) + 1

    call instantiate_patch(template_id, work_c, work_ids, coords_src, knpr_src, depth, patch)
    patch%root_element_id = root_element_id
    if (patch%n_elem > 0) then
       if (.not.allocated(patch%monitor)) allocate(patch%monitor(patch%n_elem))
       patch%monitor = monitor_value
    end if
    if (patch%n_vert <= 0) then
       call clear_patch(patch)
       return
    end if

    if (.not.patch%is_final .and. depth < max_refinement_depth) then
       if (patch%n_elem > 0 .and. allocated(patch%kvert)) then
          do i = 1, patch%n_elem
             call refine_element(root_element_id, patch%global_coor, patch%knpr, &
                  patch%kvert(:, i), monitor_value, depth + 1, used_templates, unknown_templates)
          end do
       end if
       call clear_patch(patch)
    else
       call store_final_patch(root_element_id, patch)
    end if
  end subroutine refine_element

  subroutine instantiate_patch(template_id, pattern, ivert, coords_src, knpr_src, depth, patch)
    integer, intent(in) :: template_id
    logical, intent(in) :: pattern(8)
    integer, intent(in) :: ivert(8)
    real(rk), intent(in) :: coords_src(:, :)
    integer, intent(in) :: knpr_src(:)
    integer, intent(in) :: depth
    type(element_patch_type), intent(inout) :: patch
    integer :: i, knpr_value, n_coords, n_knpr
    real(rk) :: element_coords(3, 8)

    call clear_patch(patch)
    if (template_id <= 0) return
    if (.not.allocated(meshes)) return
    if (template_id > size(meshes)) return

    associate(tmpl => meshes(template_id))
       if (tmpl%nvt <= 0 .or. .not.allocated(tmpl%coor)) return

       patch%template_id = template_id
       patch%is_final = template_is_final(template_id)
       patch%vertex_pattern = pattern
       patch%element_vertex_ids = ivert
       patch%n_vert = tmpl%nvt
       patch%n_elem = tmpl%nel
       patch%depth = depth

       allocate(patch%local_coor(3, patch%n_vert))
       patch%local_coor = tmpl%coor
       allocate(patch%global_coor(3, patch%n_vert))
       allocate(patch%kvert(8, patch%n_elem))
       if (tmpl%nel > 0 .and. allocated(tmpl%kvert)) then
          patch%kvert = tmpl%kvert
       else
          patch%kvert = 0
       end if
       allocate(patch%knpr(patch%n_vert))
       patch%knpr = 0
       if (patch%n_elem > 0) then
          allocate(patch%monitor(patch%n_elem))
          patch%monitor = 0
       end if

       n_coords = size(coords_src, 2)
       n_knpr = size(knpr_src)

       do i = 1, 8
          if (ivert(i) > 0 .and. ivert(i) <= n_coords) then
             element_coords(:, i) = coords_src(:, ivert(i))
          else
             element_coords(:, i) = 0.0_rk
          end if
          if (i <= size(patch%knpr)) then
             knpr_value = merge(1, 0, pattern(i))
             patch%knpr(i) = knpr_value
          end if
       end do

       call fill_up_element(patch, element_coords)
    end associate
  end subroutine instantiate_patch

  subroutine store_final_patch(root_element_id, patch)
    integer, intent(in) :: root_element_id
    type(element_patch_type), intent(inout) :: patch

    if (.not.allocated(element_patches)) then
      write(*, '(A)') 'Patch storage not initialized; discarding patch.'
      call clear_patch(patch)
      return
    end if

    if (root_element_id < 1 .or. root_element_id > size(element_patches)) then
       write(*, '(A,I0)') 'Invalid root element for patch storage: ', root_element_id
       call clear_patch(patch)
       return
    end if

    call append_patch_to_group(element_patches(root_element_id), patch)
    element_patch_count = element_patch_count + 1
    call clear_patch(patch)
  end subroutine store_final_patch

  subroutine fill_up_element(patch, element_coords)
    type(element_patch_type), intent(inout) :: patch
    real(rk), intent(in) :: element_coords(3, 8)

    if (patch%n_vert <= 0) return
    if (.not.allocated(patch%local_coor)) return
    if (.not.allocated(patch%global_coor)) return

    call compute_linear_coordinates(patch, element_coords, patch%global_coor)
    if (.not.apply_cylindric_transform) return

    call apply_cylindrical_blend(patch, element_coords, patch%global_coor)
  end subroutine fill_up_element

  subroutine compute_linear_coordinates(patch, element_coords, coords_out)
    type(element_patch_type), intent(in) :: patch
    real(rk), intent(in) :: element_coords(3, 8)
    real(rk), intent(inout) :: coords_out(:, :)
    real(rk), parameter :: q8 = 0.125_rk
    real(rk) :: dj(8, 3)
    real(rk) :: xi1, xi2, xi3
    real(rk) :: djac(3, 3)
    real(rk) :: xx, yy, zz
    integer :: ii
    real(rk), dimension(3, 8) :: e

    e = element_coords

    dj(1,1) = ( e(1,1)+e(1,2)+e(1,3)+e(1,4)+e(1,5)+e(1,6)+e(1,7)+e(1,8))*q8
    dj(1,2) = ( e(2,1)+e(2,2)+e(2,3)+e(2,4)+e(2,5)+e(2,6)+e(2,7)+e(2,8))*q8
    dj(1,3) = ( e(3,1)+e(3,2)+e(3,3)+e(3,4)+e(3,5)+e(3,6)+e(3,7)+e(3,8))*q8
    dj(2,1) = (-e(1,1)+e(1,2)+e(1,3)-e(1,4)-e(1,5)+e(1,6)+e(1,7)-e(1,8))*q8
    dj(2,2) = (-e(2,1)+e(2,2)+e(2,3)-e(2,4)-e(2,5)+e(2,6)+e(2,7)-e(2,8))*q8
    dj(2,3) = (-e(3,1)+e(3,2)+e(3,3)-e(3,4)-e(3,5)+e(3,6)+e(3,7)-e(3,8))*q8
    dj(3,1) = (-e(1,1)-e(1,2)+e(1,3)+e(1,4)-e(1,5)-e(1,6)+e(1,7)+e(1,8))*q8
    dj(3,2) = (-e(2,1)-e(2,2)+e(2,3)+e(2,4)-e(2,5)-e(2,6)+e(2,7)+e(2,8))*q8
    dj(3,3) = (-e(3,1)-e(3,2)+e(3,3)+e(3,4)-e(3,5)-e(3,6)+e(3,7)+e(3,8))*q8
    dj(4,1) = (-e(1,1)-e(1,2)-e(1,3)-e(1,4)+e(1,5)+e(1,6)+e(1,7)+e(1,8))*q8
    dj(4,2) = (-e(2,1)-e(2,2)-e(2,3)-e(2,4)+e(2,5)+e(2,6)+e(2,7)+e(2,8))*q8
    dj(4,3) = (-e(3,1)-e(3,2)-e(3,3)-e(3,4)+e(3,5)+e(3,6)+e(3,7)+e(3,8))*q8
    dj(5,1) = ( e(1,1)-e(1,2)+e(1,3)-e(1,4)+e(1,5)-e(1,6)+e(1,7)-e(1,8))*q8
    dj(5,2) = ( e(2,1)-e(2,2)+e(2,3)-e(2,4)+e(2,5)-e(2,6)+e(2,7)-e(2,8))*q8
    dj(5,3) = ( e(3,1)-e(3,2)+e(3,3)-e(3,4)+e(3,5)-e(3,6)+e(3,7)-e(3,8))*q8
    dj(6,1) = ( e(1,1)-e(1,2)-e(1,3)+e(1,4)-e(1,5)+e(1,6)+e(1,7)-e(1,8))*q8
    dj(6,2) = ( e(2,1)-e(2,2)-e(2,3)+e(2,4)-e(2,5)+e(2,6)+e(2,7)-e(2,8))*q8
    dj(6,3) = ( e(3,1)-e(3,2)-e(3,3)+e(3,4)-e(3,5)+e(3,6)+e(3,7)-e(3,8))*q8
    dj(7,1) = ( e(1,1)+e(1,2)-e(1,3)-e(1,4)-e(1,5)-e(1,6)+e(1,7)+e(1,8))*q8
    dj(7,2) = ( e(2,1)+e(2,2)-e(2,3)-e(2,4)-e(2,5)-e(2,6)+e(2,7)+e(2,8))*q8
    dj(7,3) = ( e(3,1)+e(3,2)-e(3,3)-e(3,4)-e(3,5)-e(3,6)+e(3,7)+e(3,8))*q8
    dj(8,1) = (-e(1,1)+e(1,2)-e(1,3)+e(1,4)+e(1,5)-e(1,6)+e(1,7)-e(1,8))*q8
    dj(8,2) = (-e(2,1)+e(2,2)-e(2,3)+e(2,4)+e(2,5)-e(2,6)+e(2,7)-e(2,8))*q8
    dj(8,3) = (-e(3,1)+e(3,2)-e(3,3)+e(3,4)+e(3,5)-e(3,6)+e(3,7)-e(3,8))*q8

    do ii = 1, min(patch%n_vert, size(coords_out, 2))
       xi1 = patch%local_coor(1, ii)
       xi2 = patch%local_coor(2, ii)
       xi3 = patch%local_coor(3, ii)

       djac(1,1) = dj(2,1) + dj(5,1)*xi2 + dj(6,1)*xi3 + dj(8,1)*xi2*xi3
       djac(1,2) = dj(3,1) + dj(5,1)*xi1 + dj(7,1)*xi3 + dj(8,1)*xi1*xi3
       djac(1,3) = dj(4,1) + dj(6,1)*xi1 + dj(7,1)*xi2 + dj(8,1)*xi1*xi2
       djac(2,1) = dj(2,2) + dj(5,2)*xi2 + dj(6,2)*xi3 + dj(8,2)*xi2*xi3
       djac(2,2) = dj(3,2) + dj(5,2)*xi1 + dj(7,2)*xi3 + dj(8,2)*xi1*xi3
       djac(2,3) = dj(4,2) + dj(6,2)*xi1 + dj(7,2)*xi2 + dj(8,2)*xi1*xi2
       djac(3,1) = dj(2,3) + dj(5,3)*xi2 + dj(6,3)*xi3 + dj(8,3)*xi2*xi3
       djac(3,2) = dj(3,3) + dj(5,3)*xi1 + dj(7,3)*xi3 + dj(8,3)*xi1*xi3
       djac(3,3) = dj(4,3) + dj(6,3)*xi1 + dj(7,3)*xi2 + dj(8,3)*xi1*xi2

       xx = dj(1,1) + djac(1,1)*xi1 + dj(3,1)*xi2 + dj(4,1)*xi3 + dj(7,1)*xi2*xi3
       yy = dj(1,2) + dj(2,2)*xi1 + djac(2,2)*xi2 + dj(4,2)*xi3 + dj(6,2)*xi1*xi3
       zz = dj(1,3) + dj(2,3)*xi1 + dj(3,3)*xi2 + djac(3,3)*xi3 + dj(5,3)*xi1*xi2

       coords_out(1, ii) = xx
       coords_out(2, ii) = yy
       coords_out(3, ii) = zz
    end do
  end subroutine compute_linear_coordinates

  subroutine apply_cylindrical_blend(patch, element_coords, coords_out)
    type(element_patch_type), intent(in) :: patch
    real(rk), intent(in) :: element_coords(3, 8)
    real(rk), intent(inout) :: coords_out(:, :)
    real(rk), allocatable :: cyl_coords(:,:)

    if (patch%n_vert <= 0) return
    if (.not.allocated(patch%local_coor)) return
    allocate(cyl_coords(3, patch%n_vert))
    call compute_cylindrical_coordinates(patch, element_coords, cyl_coords)
    call blend_coordinate_sets(coords_out, cyl_coords)
    deallocate(cyl_coords)
  end subroutine apply_cylindrical_blend
  subroutine mark_elements_by_threshold(mesh, level, threshold_value)
    type(mesh_type), intent(inout) :: mesh
    integer, intent(in) :: level
    real(rk), intent(in) :: threshold_value
    integer :: elem, max_elem

    if (.not.allocated(mesh%monitor)) return
    if (.not.allocated(mesh%monitor_value)) return
    max_elem = min(size(mesh%monitor), size(mesh%monitor_value))

    do elem = 1, max_elem
       if (mesh%monitor_value(elem) > threshold_value) then
          if (mesh%monitor(elem) < level) mesh%monitor(elem) = level
       end if
    end do
  end subroutine mark_elements_by_threshold

  subroutine enforce_refinement_levels(mesh, level)
    type(mesh_type), intent(inout) :: mesh
    integer, intent(in) :: level
    integer :: elem, neighbor
    logical :: changed

    if (.not.allocated(mesh%monitor)) return
    if (.not.allocated(mesh%kelementspan)) return

    do
       changed = .false.
       do elem = 1, min(size(mesh%monitor), size(mesh%kelementspan))
          if (mesh%monitor(elem) /= level) cycle
          if (.not.allocated(mesh%kelementspan(elem)%list)) cycle
          do neighbor = 1, size(mesh%kelementspan(elem)%list)
             if (mesh%kelementspan(elem)%list(neighbor) < 1) cycle
             if (mesh%kelementspan(elem)%list(neighbor) > size(mesh%monitor)) cycle
             if (mesh%monitor(mesh%kelementspan(elem)%list(neighbor)) < level - 1) then
                mesh%monitor(mesh%kelementspan(elem)%list(neighbor)) = level - 1
                changed = .true.
             end if
          end do
       end do
       if (.not.changed) exit
    end do
  end subroutine enforce_refinement_levels

  subroutine apply_inflow_refinement_boost(mesh, setup_path, changed)
    type(mesh_type), intent(inout) :: mesh
    character(len=*), intent(in) :: setup_path
    logical, intent(out) :: changed

    type(MeshConfig) :: mesh_config
    type(ProcessParameters) :: process_params
    type(FaceList) :: boundary_faces
    type(BoxBoundaryClassification) :: box_boundary
    type(HollowCylinderBoundaryClassification) :: cyl_boundary
    real(c_double), allocatable :: coords(:, :)
    integer(c_int), allocatable :: kvert(:, :)
    integer(c_int), allocatable :: knpr(:)
    logical, allocatable :: inflow_mask(:)
    integer :: nel, nvt, elem, boosted
    logical :: file_exists, is_box_mesh, is_cyl_mesh
    character(len=64) :: type_token

    changed = .false.
    inquire(file=trim(setup_path), exist=file_exists)
    if (.not.file_exists) return
    if (mesh%nel <= 0 .or. mesh%nvt <= 0) return
    if (.not.allocated(mesh%monitor)) return

    call initialize_mesh_config(mesh_config)
    call load_mesh_config(trim(setup_path), mesh_config)
    if (.not.mesh_config%loaded) return

    call initialize_process_parameters(process_params)
    call load_process_parameters(trim(setup_path), process_params)
    if (process_params%nOfInflows <= 0) return
    if (.not.allocated(process_params%inflows)) return

    type_token = lowercase_string(trim(mesh_config%mesh_type))
    is_box_mesh = (type_token == 'box')
    is_cyl_mesh = (type_token == 'hollowcylinder' .or. type_token == 'fullcylinder')
    if (.not.(is_box_mesh .or. is_cyl_mesh)) then
       write(*, '(A)') 'Unsupported HexMesher type in setup.e3d: ' // trim(mesh_config%mesh_type)
       stop 1
    end if

    nvt = mesh%nvt
    nel = mesh%nel
    allocate(coords(3, nvt))
    coords = real(mesh%coor, kind=c_double)
    allocate(kvert(8, nel))
    kvert = int(mesh%kvert, kind=c_int)
    allocate(knpr(nvt))
    knpr = 0_c_int
    allocate(inflow_mask(nel))
    inflow_mask = .false.

    call recompute_knpr_from_connectivity(kvert, knpr, boundary_faces)

    if (is_box_mesh) then
       call classify_box_boundaries(mesh_config, process_params, coords, kvert, knpr, boundary_faces, box_boundary)
       if (allocated(box_boundary%inflow_groups)) then
          call mark_inflow_elements(box_boundary%inflow_groups, inflow_mask)
       end if
    else
        call classify_hollowcylinder_boundaries(mesh_config, process_params, coords, kvert, knpr, boundary_faces, &
             cyl_boundary)
        if (allocated(cyl_boundary%inflow_groups)) then
           call mark_inflow_elements(cyl_boundary%inflow_groups, inflow_mask)
        end if
    end if

    if (.not.any(inflow_mask)) then
       call cleanup_inflow_buffers()
       return
    end if

    if (.not.allocated(mesh%kadj)) call build_face_adjacency(mesh)

    boosted = 0
    do elem = 1, min(nel, size(mesh%monitor))
       if (.not.inflow_mask(elem)) cycle
       if (.not.element_is_boundary(mesh, elem)) cycle
       if (mesh%monitor(elem) >= refinement_depth) cycle
       mesh%monitor(elem) = mesh%monitor(elem) + 1
       boosted = boosted + 1
    end do

    if (boosted > 0) then
       changed = .true.
       write(*, '(A,I0)') 'Inflow refinement boost applied to elements: ', boosted
    else
       write(*, '(A)') 'No boundary elements matched inflow definitions for refinement boost.'
    end if

    call cleanup_inflow_buffers()

  contains

    subroutine cleanup_inflow_buffers()
      if (allocated(coords)) deallocate(coords)
      if (allocated(kvert))  deallocate(kvert)
      if (allocated(knpr))   deallocate(knpr)
      if (allocated(inflow_mask)) deallocate(inflow_mask)
    end subroutine cleanup_inflow_buffers

  end subroutine apply_inflow_refinement_boost

  subroutine mark_inflow_elements(groups, mask)
    type(InflowBoundaryGroup), intent(in) :: groups(:)
    logical, intent(inout) :: mask(:)
    integer :: inflow_idx, face_idx, elem

    if (size(mask) <= 0) return

    do inflow_idx = 1, size(groups)
       if (groups(inflow_idx)%faces%count <= 0) cycle
       do face_idx = 1, groups(inflow_idx)%faces%count
          elem = groups(inflow_idx)%faces%elements(face_idx)
          if (elem >= 1 .and. elem <= size(mask)) mask(elem) = .true.
       end do
    end do
  end subroutine mark_inflow_elements

  logical function element_is_boundary(mesh, elem)
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: elem
    integer :: face, max_faces

    element_is_boundary = .false.
    if (.not.allocated(mesh%kadj)) return
    if (elem < 1 .or. elem > size(mesh%kadj, 2)) return

    max_faces = min(size(mesh%kadj, 1), 6)
    do face = 1, max_faces
       if (mesh%kadj(face, elem) <= 0) then
          element_is_boundary = .true.
          return
       end if
    end do
  end function element_is_boundary

  pure function lowercase_string(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: idx, code

    out = text
    do idx = 1, len(text)
       code = iachar(text(idx:idx))
       if (code >= iachar('A') .and. code <= iachar('Z')) then
          out(idx:idx) = achar(code + 32)
       else
          out(idx:idx) = text(idx:idx)
       end if
    end do
  end function lowercase_string

  subroutine report_refinement_distribution(mesh, label)
    type(mesh_type), intent(in) :: mesh
    character(len=*), intent(in) :: label
    integer :: elem, lvl
    integer, allocatable :: counts(:)

    if (.not.allocated(mesh%monitor)) then
       write(*, '(A)') trim(label) // ': monitor data not available.'
       return
    end if
    allocate(counts(0:refinement_depth))
    counts = 0

    do elem = 1, min(size(mesh%monitor), mesh%nel)
       lvl = mesh%monitor(elem)
       if (lvl < 0) lvl = 0
       if (lvl > refinement_depth) lvl = refinement_depth
       counts(lvl) = counts(lvl) + 1
    end do

    write(*, '(A)', advance='no') trim(label) // ' ['
    do lvl = 0, refinement_depth
       if (lvl > 0) write(*, '(A)', advance='no') ', '
       write(*, '(I0)', advance='no') counts(lvl)
    end do
    write(*, '(A)') ']'

    deallocate(counts)
  end subroutine report_refinement_distribution

  subroutine compute_cylindrical_coordinates(patch, element_coords, coords_out)
    type(element_patch_type), intent(in) :: patch
    real(rk), intent(in) :: element_coords(3, 8)
    real(rk), intent(inout) :: coords_out(:, :)
    real(rk), parameter :: q8 = 0.125_rk
    real(rk) :: dj(8, 3)
    real(rk) :: xi1, xi2, xi3
    real(rk) :: djac(3, 3)
    real(rk) :: rr, tt, zz, xx, yy
    real(rk) :: ecyl(3, 8)
    real(rk) :: dth, th0
    real(rk), parameter :: pi = 4.0_rk * atan(1.0_rk)
    real(rk), parameter :: twopi = 2.0_rk * pi
    integer :: ii

    if (patch%n_vert <= 0) return
    if (.not.allocated(patch%local_coor)) return

    do ii = 1, 8
       ecyl(1, ii) = sqrt(element_coords(1, ii)**2 + element_coords(2, ii)**2)
       ecyl(2, ii) = atan2(element_coords(2, ii), element_coords(1, ii))
       ecyl(3, ii) = element_coords(3, ii)
    end do

    th0 = ecyl(2, 1)
    do ii = 2, 8
       dth = ecyl(2, ii) - th0
       if (dth > pi)  ecyl(2, ii) = ecyl(2, ii) - twopi
       if (dth < -pi) ecyl(2, ii) = ecyl(2, ii) + twopi
    end do

    dj(1,1)=( ecyl(1,1)+ecyl(1,2)+ecyl(1,3)+ecyl(1,4)+ecyl(1,5)+ecyl(1,6)+ecyl(1,7)+ecyl(1,8))*q8
    dj(1,2)=( ecyl(2,1)+ecyl(2,2)+ecyl(2,3)+ecyl(2,4)+ecyl(2,5)+ecyl(2,6)+ecyl(2,7)+ecyl(2,8))*q8
    dj(1,3)=( ecyl(3,1)+ecyl(3,2)+ecyl(3,3)+ecyl(3,4)+ecyl(3,5)+ecyl(3,6)+ecyl(3,7)+ecyl(3,8))*q8
    dj(2,1)=(-ecyl(1,1)+ecyl(1,2)+ecyl(1,3)-ecyl(1,4)-ecyl(1,5)+ecyl(1,6)+ecyl(1,7)-ecyl(1,8))*q8
    dj(2,2)=(-ecyl(2,1)+ecyl(2,2)+ecyl(2,3)-ecyl(2,4)-ecyl(2,5)+ecyl(2,6)+ecyl(2,7)-ecyl(2,8))*q8
    dj(2,3)=(-ecyl(3,1)+ecyl(3,2)+ecyl(3,3)-ecyl(3,4)-ecyl(3,5)+ecyl(3,6)+ecyl(3,7)-ecyl(3,8))*q8
    dj(3,1)=(-ecyl(1,1)-ecyl(1,2)+ecyl(1,3)+ecyl(1,4)-ecyl(1,5)-ecyl(1,6)+ecyl(1,7)+ecyl(1,8))*q8
    dj(3,2)=(-ecyl(2,1)-ecyl(2,2)+ecyl(2,3)+ecyl(2,4)-ecyl(2,5)-ecyl(2,6)+ecyl(2,7)+ecyl(2,8))*q8
    dj(3,3)=(-ecyl(3,1)-ecyl(3,2)+ecyl(3,3)+ecyl(3,4)-ecyl(3,5)-ecyl(3,6)+ecyl(3,7)+ecyl(3,8))*q8
    dj(4,1)=(-ecyl(1,1)-ecyl(1,2)-ecyl(1,3)-ecyl(1,4)+ecyl(1,5)+ecyl(1,6)+ecyl(1,7)+ecyl(1,8))*q8
    dj(4,2)=(-ecyl(2,1)-ecyl(2,2)-ecyl(2,3)-ecyl(2,4)+ecyl(2,5)+ecyl(2,6)+ecyl(2,7)+ecyl(2,8))*q8
    dj(4,3)=(-ecyl(3,1)-ecyl(3,2)-ecyl(3,3)-ecyl(3,4)+ecyl(3,5)+ecyl(3,6)+ecyl(3,7)+ecyl(3,8))*q8
    dj(5,1)=( ecyl(1,1)-ecyl(1,2)+ecyl(1,3)-ecyl(1,4)+ecyl(1,5)-ecyl(1,6)+ecyl(1,7)-ecyl(1,8))*q8
    dj(5,2)=( ecyl(2,1)-ecyl(2,2)+ecyl(2,3)-ecyl(2,4)+ecyl(2,5)-ecyl(2,6)+ecyl(2,7)-ecyl(2,8))*q8
    dj(5,3)=( ecyl(3,1)-ecyl(3,2)+ecyl(3,3)-ecyl(3,4)+ecyl(3,5)-ecyl(3,6)+ecyl(3,7)-ecyl(3,8))*q8
    dj(6,1)=( ecyl(1,1)-ecyl(1,2)-ecyl(1,3)+ecyl(1,4)-ecyl(1,5)+ecyl(1,6)+ecyl(1,7)-ecyl(1,8))*q8
    dj(6,2)=( ecyl(2,1)-ecyl(2,2)-ecyl(2,3)+ecyl(2,4)-ecyl(2,5)+ecyl(2,6)+ecyl(2,7)-ecyl(2,8))*q8
    dj(6,3)=( ecyl(3,1)-ecyl(3,2)-ecyl(3,3)+ecyl(3,4)-ecyl(3,5)+ecyl(3,6)+ecyl(3,7)-ecyl(3,8))*q8
    dj(7,1)=( ecyl(1,1)+ecyl(1,2)-ecyl(1,3)-ecyl(1,4)-ecyl(1,5)-ecyl(1,6)+ecyl(1,7)+ecyl(1,8))*q8
    dj(7,2)=( ecyl(2,1)+ecyl(2,2)-ecyl(2,3)-ecyl(2,4)-ecyl(2,5)-ecyl(2,6)+ecyl(2,7)+ecyl(2,8))*q8
    dj(7,3)=( ecyl(3,1)+ecyl(3,2)-ecyl(3,3)-ecyl(3,4)-ecyl(3,5)-ecyl(3,6)+ecyl(3,7)+ecyl(3,8))*q8
    dj(8,1)=(-ecyl(1,1)+ecyl(1,2)-ecyl(1,3)+ecyl(1,4)+ecyl(1,5)-ecyl(1,6)+ecyl(1,7)-ecyl(1,8))*q8
    dj(8,2)=(-ecyl(2,1)+ecyl(2,2)-ecyl(2,3)+ecyl(2,4)+ecyl(2,5)-ecyl(2,6)+ecyl(2,7)-ecyl(2,8))*q8
    dj(8,3)=(-ecyl(3,1)+ecyl(3,2)-ecyl(3,3)+ecyl(3,4)+ecyl(3,5)-ecyl(3,6)+ecyl(3,7)-ecyl(3,8))*q8

    do ii = 1, min(patch%n_vert, size(coords_out, 2))
       xi1 = patch%local_coor(1, ii)
       xi2 = patch%local_coor(2, ii)
       xi3 = patch%local_coor(3, ii)

       djac(1,1)=dj(2,1) + dj(5,1)*xi2 + dj(6,1)*xi3 + dj(8,1)*xi2*xi3
       djac(1,2)=dj(3,1) + dj(5,1)*xi1 + dj(7,1)*xi3 + dj(8,1)*xi1*xi3
       djac(1,3)=dj(4,1) + dj(6,1)*xi1 + dj(7,1)*xi2 + dj(8,1)*xi1*xi2
       djac(2,1)=dj(2,2) + dj(5,2)*xi2 + dj(6,2)*xi3 + dj(8,2)*xi2*xi3
       djac(2,2)=dj(3,2) + dj(5,2)*xi1 + dj(7,2)*xi3 + dj(8,2)*xi1*xi3
       djac(2,3)=dj(4,2) + dj(6,2)*xi1 + dj(7,2)*xi2 + dj(8,2)*xi1*xi2
       djac(3,1)=dj(2,3) + dj(5,3)*xi2 + dj(6,3)*xi3 + dj(8,3)*xi2*xi3
       djac(3,2)=dj(3,3) + dj(5,3)*xi1 + dj(7,3)*xi3 + dj(8,3)*xi1*xi3
       djac(3,3)=dj(4,3) + dj(6,3)*xi1 + dj(7,3)*xi2 + dj(8,3)*xi1*xi2

       rr = dj(1,1) + djac(1,1)*xi1 + dj(3,1)*xi2 + dj(4,1)*xi3 + dj(7,1)*xi2*xi3
       tt = dj(1,2) + dj(2,2)*xi1 + djac(2,2)*xi2 + dj(4,2)*xi3 + dj(6,2)*xi1*xi3
       zz = dj(1,3) + dj(2,3)*xi1 + dj(3,3)*xi2 + djac(3,3)*xi3 + dj(5,3)*xi1*xi2

       xx = rr * cos(tt)
       yy = rr * sin(tt)

        coords_out(1, ii) = xx
        coords_out(2, ii) = yy
        coords_out(3, ii) = zz
    end do
  end subroutine compute_cylindrical_coordinates

  subroutine blend_coordinate_sets(box_coords, cyl_coords)
    real(rk), intent(inout) :: box_coords(:, :)
    real(rk), intent(in) :: cyl_coords(:, :)
    real(rk) :: outer_radius, diameter, low_thresh, high_thresh, radius, alpha
    real(rk) :: pct_min, pct_max
    integer :: ii, max_nodes

    max_nodes = min(size(box_coords, 2), size(cyl_coords, 2))
    if (max_nodes <= 0) return

    outer_radius = cylindrical_outer_radius
    if (outer_radius <= 0.0_rk) then
       outer_radius = 0.0_rk
       do ii = 1, max_nodes
          radius = sqrt(box_coords(1, ii)**2 + box_coords(2, ii)**2)
          if (radius > outer_radius) outer_radius = radius
       end do
    end if

    diameter = 2.0_rk * outer_radius
    if (diameter <= 0.0_rk) then
       do ii = 1, max_nodes
          box_coords(:, ii) = cyl_coords(:, ii)
       end do
       return
    end if

    pct_min = max(0.0_rk, real(get_fullcyl_min_radius_percentage(), rk))
    pct_max = max(0.0_rk, real(get_fullcyl_max_radius_percentage(), rk))
    if (pct_max <= pct_min) pct_max = pct_min + 1.0_rk

    low_thresh = 0.01_rk * pct_min * diameter
    high_thresh = 0.01_rk * pct_max * diameter

    do ii = 1, max_nodes
       radius = sqrt(box_coords(1, ii)**2 + box_coords(2, ii)**2)
       if (radius <= low_thresh) then
          alpha = 0.0_rk
       else if (radius >= high_thresh) then
          alpha = 1.0_rk
       else
          alpha = (radius - low_thresh) / (high_thresh - low_thresh)
       end if
       box_coords(:, ii) = alpha * cyl_coords(:, ii) + (1.0_rk - alpha) * box_coords(:, ii)
    end do
  end subroutine blend_coordinate_sets

  subroutine write_template_usage(usage)
    integer, intent(in) :: usage(:)
    integer :: i
    character(len=16) :: buf
    character(len=3) :: status

    write(*, '(A)', advance='no') 'Templates '
    do i = 1, size(usage)
       write(buf, '(I0)') i
       write(*, '(A)', advance='no') '[' // pad_field(buf) // ']'
    end do
    write(*, '(A)') ''

    write(*, '(A)', advance='no') 'Status    '
    do i = 1, size(usage)
       status = merge('FIN', 'INT', template_is_final(i))
       write(*, '(A)', advance='no') '[' // pad_field(status) // ']'
    end do
    write(*, '(A)') ''

    write(*, '(A)', advance='no') 'Usage     '
    do i = 1, size(usage)
       write(buf, '(I0)') usage(i)
       write(*, '(A)', advance='no') '[' // pad_field(buf) // ']'
    end do
    write(*, '(A)') ''
  contains
    pure function pad_field(text) result(padded)
      character(len=*), intent(in) :: text
      character(len=8) :: padded
      integer :: len_text, start_idx

      padded = '        '
      len_text = len_trim(text)
      if (len_text <= 0) return
      if (len_text >= len(padded)) then
         padded = text(len_text - len(padded) + 1:len_text)
      else
         start_idx = len(padded) - len_text + 1
         padded(start_idx:) = text(1:len_text)
      end if
    end function pad_field
  end subroutine write_template_usage

end module def_mod
