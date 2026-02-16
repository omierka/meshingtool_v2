module var_mod
  implicit none

  integer, parameter :: rk = kind(1.0d0)
  integer, parameter :: default_refinement_depth = 2
  integer :: level
  real(rk) :: Monitor_threshold = 1.5_rk
  real(rk) :: scaling_factor_TRI_output = 0.1_rk


  integer, parameter :: mesh_file_count = 23
  character(len=32), parameter :: mesh_files(mesh_file_count) = (/ &
         'VERTEX/FIN/00000000.tri', &
         'VERTEX/FIN/10000000.tri', &
         'VERTEX/FIN/10000010.tri', &
         'VERTEX/FIN/10100000.tri', &
         'VERTEX/FIN/10100100.tri', &
         'VERTEX/FIN/10100101.tri', &
         'VERTEX/FIN/10101000.tri', &
         'VERTEX/FIN/11000000.tri', &
         'VERTEX/FIN/11110000.tri', &
         'VERTEX/FIN/11111111.tri', &
         'VERTEX/INT/10101010.tri', &
         'VERTEX/INT/10110100.tri', &
         'VERTEX/INT/10111000.tri', &
         'VERTEX/INT/10111100.tri', &
         'VERTEX/INT/10111110.tri', &
         'VERTEX/INT/11011000.tri', &
         'VERTEX/INT/11011010.tri', &
         'VERTEX/INT/11100000.tri', &
         'VERTEX/INT/11101000.tri', &
         'VERTEX/INT/11111000.tri', &
         'VERTEX/INT/11111010.tri', &
         'VERTEX/INT/11111100.tri', &
         'VERTEX/INT/11111110.tri' /)
  logical :: template_is_final(mesh_file_count)
  data template_is_final / mesh_file_count*.false. /

  logical, parameter :: templates(8, mesh_file_count) = reshape([ &
       .false., .false., .false., .false., .false., .false., .false., .false., &
       .true.,  .false., .false., .false., .false., .false., .false., .false., &
       .true.,  .true.,  .false., .false., .false., .false., .false., .false., &
       .true.,  .false., .true.,  .false., .false., .false., .false., .false., &
       .true.,  .false., .false., .false., .false., .false., .true.,  .false., &
       .true.,  .true.,  .true.,  .false., .false., .false., .false., .false., &
       .true.,  .false., .true.,  .false., .true.,  .false., .false., .false., &
       .true.,  .false., .true.,  .false., .false., .true.,  .false., .false., &
       .true.,  .true.,  .true.,  .true.,  .false., .false., .false., .false., &
       .true.,  .true.,  .false., .true.,  .true.,  .false., .false., .false., &
       .true.,  .false., .true.,  .true.,  .true.,  .false., .false., .false., &
       .true.,  .false., .true.,  .true.,  .false., .true.,  .false., .false., &
       .true.,  .false., .true.,  .false., .true.,  .false., .true.,  .false., &
       .true.,  .false., .true.,  .false., .false., .true.,  .false., .true.,  &
       .true.,  .true.,  .true.,  .true.,  .true.,  .false., .false., .false., &
       .true.,  .false., .true.,  .true.,  .true.,  .true.,  .false., .false., &
       .true.,  .true.,  .false., .true.,  .true.,  .false., .true.,  .false., &
       .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .false., .false., &
       .true.,  .true.,  .true.,  .true.,  .true.,  .false., .true.,  .false., &
       .true.,  .false., .true.,  .true.,  .true.,  .true.,  .true.,  .false., &
       .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .false., &
       .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  &
       .true.,  .true.,  .true.,  .false., .true.,  .false., .false., .false.  ], &
       shape=[8, mesh_file_count])

  integer :: RROT(8,8), R2ROT(8), R3ROT(8)
  data RROT/1,2,3,4,5,6,7,8,&
         2,3,4,1,6,7,8,5,&
         3,4,1,2,7,8,5,6,&
         4,1,2,3,8,5,6,7,&
         5,8,7,6,1,4,3,2,&
         6,5,8,7,2,1,4,3,&
         7,6,5,8,3,2,1,4,&
         8,7,6,5,4,3,2,1/
  data R2ROT/1,5,6,2,4,8,7,3/
  data R3ROT/1,4,8,5,2,3,7,6/

  type :: element_span_type
     integer, allocatable :: list(:)
  end type element_span_type

  type :: mesh_type
     integer :: nel  = 0
     integer :: nvt  = 0
     integer :: nbct = 0
     integer :: nve  = 0
     integer :: nee  = 0
     integer :: nae  = 0
     integer :: template_id = 0
     real(rk), allocatable    :: coor(:,:)  !! Cartesian coordinates, size (3,nvt)
     integer, allocatable     :: kvert(:,:) !! Vertex connectivities, size (8,nel)
     integer, allocatable     :: knpr(:)    !! Nodal properties, size (nvt)
     integer, allocatable     :: monitor(:) !! Element tags, size (nel)
     real(rk), allocatable    :: monitor_value(:) !! Raw monitor values
     integer, allocatable     :: kadj(:,:)  !! Element adjacency, size (6,nel)
     type(element_span_type), allocatable :: kelementspan(:)
  end type mesh_type

  type :: element_patch_type
     integer :: template_id = 0
     logical :: is_final = .false.
     logical :: vertex_pattern(8) = .false.
     integer :: element_vertex_ids(8) = 0
     integer :: n_vert = 0
     integer :: n_elem = 0
     integer :: depth = 0
     integer :: root_element_id = 0
     real(rk), allocatable :: local_coor(:,:)
     real(rk), allocatable :: global_coor(:,:)
     integer, allocatable :: kvert(:,:)
     integer, allocatable :: knpr(:)
     integer, allocatable :: monitor(:)
  end type element_patch_type

  type :: element_patch_group
     type(element_patch_type), allocatable :: patchlist(:)
     integer :: count = 0
  end type element_patch_group

  type(mesh_type), allocatable :: meshes(:)
  type(mesh_type), allocatable, target :: hex_mesh(:)
  type(mesh_type), pointer :: target_mesh => null()
  type(mesh_type), pointer :: refined_clean_mesh => null()
  integer :: refinement_depth = 0
  logical :: reproducibility = .true.
  logical :: apply_cylindric_transform = .false.
  type(element_patch_group), allocatable, target :: element_patches(:)
  type(element_patch_group), allocatable, target :: clean_element_patches(:)
  integer :: element_patch_count = 0

contains

  subroutine clear_mesh(mesh)
    type(mesh_type), intent(inout) :: mesh

    if (allocated(mesh%coor))  deallocate(mesh%coor)
    if (allocated(mesh%kvert)) deallocate(mesh%kvert)
    if (allocated(mesh%knpr))  deallocate(mesh%knpr)
    if (allocated(mesh%kadj))  deallocate(mesh%kadj)
    if (allocated(mesh%monitor)) deallocate(mesh%monitor)
    if (allocated(mesh%monitor_value)) deallocate(mesh%monitor_value)
    if (allocated(mesh%kelementspan)) then
       call release_element_span(mesh%kelementspan)
       deallocate(mesh%kelementspan)
    end if

    mesh%nel  = 0
    mesh%nvt  = 0
    mesh%nbct = 0
    mesh%nve  = 0
    mesh%nee  = 0
    mesh%nae  = 0
    mesh%template_id = 0
  end subroutine clear_mesh

  subroutine clear_patch(patch)
    type(element_patch_type), intent(inout) :: patch

    patch%template_id = 0
    patch%is_final = .false.
    patch%vertex_pattern = .false.
    patch%element_vertex_ids = 0
     patch%n_vert = 0
     patch%n_elem = 0
     patch%depth = 0
    patch%root_element_id = 0
     if (allocated(patch%local_coor))  deallocate(patch%local_coor)
     if (allocated(patch%global_coor)) deallocate(patch%global_coor)
     if (allocated(patch%kvert))       deallocate(patch%kvert)
    if (allocated(patch%knpr))        deallocate(patch%knpr)
    if (allocated(patch%monitor))     deallocate(patch%monitor)
  end subroutine clear_patch

  subroutine release_element_span(span_array)
    type(element_span_type), allocatable, intent(inout) :: span_array(:)
    integer :: i

    if (.not.allocated(span_array)) return
    do i = 1, size(span_array)
       if (allocated(span_array(i)%list)) deallocate(span_array(i)%list)
    end do
  end subroutine release_element_span

  subroutine release_element_patches()
    integer :: i

    element_patch_count = 0
    if (allocated(element_patches)) then
       do i = 1, size(element_patches)
          call clear_patch_group(element_patches(i))
       end do
       deallocate(element_patches)
    end if
    if (allocated(clean_element_patches)) then
       do i = 1, size(clean_element_patches)
          call clear_patch_group(clean_element_patches(i))
       end do
       deallocate(clean_element_patches)
    end if
  end subroutine release_element_patches

  subroutine initialize_mesh_levels(depth)
    integer, intent(in) :: depth
    integer :: i

    call release_mesh_levels()
    allocate(hex_mesh(0:depth))
    do i = 0, depth
       call clear_mesh(hex_mesh(i))
    end do
    refinement_depth = depth
    target_mesh => hex_mesh(0)
    if (depth >= 1) then
       refined_clean_mesh => hex_mesh(1)
    else
       nullify(refined_clean_mesh)
    end if
  end subroutine initialize_mesh_levels

  subroutine release_mesh_levels()
    integer :: i

    if (.not.allocated(hex_mesh)) then
       nullify(target_mesh)
       nullify(refined_clean_mesh)
       refinement_depth = 0
       return
    end if
    do i = lbound(hex_mesh, 1), ubound(hex_mesh, 1)
       call clear_mesh(hex_mesh(i))
    end do
    deallocate(hex_mesh)
    nullify(target_mesh)
    nullify(refined_clean_mesh)
    refinement_depth = 0
  end subroutine release_mesh_levels

  subroutine bind_refined_mesh(level)
    integer, intent(in) :: level

    if (.not.allocated(hex_mesh)) then
      nullify(refined_clean_mesh)
      return
    end if
    if (level < lbound(hex_mesh, 1) .or. level > ubound(hex_mesh, 1)) then
      nullify(refined_clean_mesh)
    else
      refined_clean_mesh => hex_mesh(level)
    end if
  end subroutine bind_refined_mesh

  subroutine clear_patch_group(group)
    type(element_patch_group), intent(inout) :: group
    integer :: i

    if (allocated(group%patchlist)) then
       do i = 1, size(group%patchlist)
          call clear_patch(group%patchlist(i))
       end do
       deallocate(group%patchlist)
    end if
    group%count = 0
  end subroutine clear_patch_group

  subroutine append_patch_to_group(group, patch)
    type(element_patch_group), intent(inout) :: group
    type(element_patch_type), intent(in) :: patch
    type(element_patch_type), allocatable :: newlist(:)
    integer :: new_count

    new_count = group%count + 1
    if (.not.allocated(group%patchlist)) then
       allocate(group%patchlist(new_count))
    else
       allocate(newlist(new_count))
       newlist(1:group%count) = group%patchlist
       call move_alloc(newlist, group%patchlist)
    end if
    group%patchlist(new_count) = patch
    group%count = new_count
  end subroutine append_patch_to_group

  pure function extract_template_code(path) result(code)
    character(len=*), intent(in) :: path
    character(len=8) :: code
    integer :: pos, digit_count
    character(len=1) :: ch

    code = '00000000'
    digit_count = 0
    do pos = len_trim(path), 1, -1
       ch = path(pos:pos)
       if (ch == '0' .or. ch == '1') then
          digit_count = digit_count + 1
          if (digit_count <= 8) then
             code(9 - digit_count:9 - digit_count) = ch
          end if
          if (digit_count == 8) exit
       end if
    end do
  end function extract_template_code

  pure subroutine code_to_pattern(code, pattern)
    character(len=8), intent(in) :: code
    logical, intent(out) :: pattern(8)
    integer :: j

    do j = 1, 8
       pattern(j) = (code(j:j) == '1')
    end do
  end subroutine code_to_pattern

  integer function match_template_pattern(pattern) result(idx)
    logical, intent(in) :: pattern(8)
    logical :: work_c(8)
    integer :: dummy_ids(8)
    integer :: i

    idx = 0
    do i = 1, mesh_file_count
       if (all(pattern .eqv. templates(:, i))) then
          idx = i
          return
       end if
    end do

    work_c = pattern
    dummy_ids = [(i, i = 1, 8)]
    call rotate_patch(work_c, dummy_ids)

    do i = 1, mesh_file_count
       if (all(work_c .eqv. templates(:, i))) then
          idx = i
          return
       end if
    end do
  end function match_template_pattern

  integer function template_index_from_filename(path) result(idx)
    character(len=*), intent(in) :: path
    character(len=8) :: code
    logical :: pattern(8)

    code = extract_template_code(path)
    call code_to_pattern(code, pattern)
    idx = match_template_pattern(pattern)
  end function template_index_from_filename

  pure function pattern_to_string(pattern) result(code)
    logical, intent(in) :: pattern(8)
    character(len=8) :: code
    integer :: j

    do j = 1, 8
       code(j:j) = merge('1', '0', pattern(j))
    end do
  end function pattern_to_string

  subroutine rotate_patch(c, id)
    logical, intent(inout) :: c(8)
    integer, intent(inout) :: id(8)
    logical :: r1(8), r2(8), r3(8), p(8)
    integer :: idr(8), idr1(8), idr2(8), idr3(8)
    integer :: i, j, isum, jsum, jsum_init

    jsum_init = 64 * 8
    jsum = jsum_init

    do i = 1, 8
       if (c(i)) then
          r1 = c(rrot(:, i))
          idr1 = id(rrot(:, i))
          isum = 0
          do j = 1, 8
             if (r1(j)) isum = isum + j * j
          end do
          if (isum < jsum) then
             p = r1
             idr = idr1
             jsum = isum
          end if

          r2 = r1(r2rot)
          idr2 = idr1(r2rot)
          isum = 0
          do j = 1, 8
             if (r2(j)) isum = isum + j * j
          end do
          if (isum < jsum) then
             p = r2
             idr = idr2
             jsum = isum
          end if

          r3 = r1(r3rot)
          idr3 = idr1(r3rot)
          isum = 0
          do j = 1, 8
             if (r3(j)) isum = isum + j * j
          end do
          if (isum < jsum) then
             p = r3
             idr = idr3
             jsum = isum
          end if
       end if
    end do

    if (jsum < jsum_init) then
       id = idr
       c = p
    end if
  end subroutine rotate_patch

  subroutine determine_template(c, it, bfinal)
    logical, intent(in)  :: c(8)
    integer, intent(out) :: it
    logical, intent(out), optional :: bfinal
    logical :: r(8), bsame
    integer :: i

    it = 0
    if (present(bfinal)) bfinal = .false.
    do i = 1, size(templates, 2)
       r = templates(:, i)
       bsame = .true.
       call compare_patch(c, r, bsame)
       if (.not.bsame) then
          it = i
          if (present(bfinal)) bfinal = template_is_final(i)
          return
       end if
    end do
  end subroutine determine_template

  subroutine compare_patch(c, r, b)
    logical, intent(in)    :: c(8), r(8)
    logical, intent(inout) :: b
    integer :: i, match_count

    match_count = 0
    do i = 1, 8
       if ((c(i) .and. r(i)) .or. ((.not.c(i)) .and. (.not.r(i)))) match_count = match_count + 1
    end do

    if (match_count == 8) b = .false.
  end subroutine compare_patch

end module var_mod
