module inout_mod
  use var_mod, only: mesh_type, meshes, clear_mesh, rk, mesh_file_count, mesh_files, target_mesh, &
                     extract_template_code, code_to_pattern, match_template_pattern, pattern_to_string, &
                     rotate_patch, templates, template_is_final, element_patch_group, apply_cylindric_transform, &
                     scaling_factor_TRI_output, cylindrical_outer_radius
  implicit none

contains

  subroutine load_all_meshes()
    integer :: i, j, idx
    character(len=8) :: code, normalized_code, canonical_code
    logical :: pattern(8), rotated_pattern(8)
    integer :: dummy_ids(8)
    logical :: has_failure, is_final
    logical, allocatable :: template_loaded(:)
    character(len=32) :: buf

    if (allocated(meshes)) call release_all_meshes()
    allocate(meshes(mesh_file_count))
    do i = 1, mesh_file_count
       call clear_mesh(meshes(i))
    end do
    template_is_final = .false.
    allocate(template_loaded(mesh_file_count))
    template_loaded = .false.

    has_failure = .false.
    do i = 1, mesh_file_count
       code = extract_template_code(mesh_files(i))
       call code_to_pattern(code, pattern)
       idx = match_template_pattern(pattern)
       if (idx <= 0 .or. idx > mesh_file_count) then
          rotated_pattern = pattern
          dummy_ids = [(j, j = 1, 8)]
          call rotate_patch(rotated_pattern, dummy_ids)
          normalized_code = pattern_to_string(rotated_pattern)
          write(*, '(A)') 'Unable to map template file to canonical template: ' // trim(mesh_files(i))
          write(*, '(A,A)') ' raw pattern: ', code
          write(*, '(A,A)') ' rotated pattern: ', normalized_code
          has_failure = .true.
          cycle
       end if
       rotated_pattern = pattern
       dummy_ids = [(j, j = 1, 8)]
       call rotate_patch(rotated_pattern, dummy_ids)
       normalized_code = pattern_to_string(rotated_pattern)
       canonical_code = pattern_to_string(templates(:, idx))
      call clear_mesh(meshes(idx))
       call read_mesh_file(mesh_files(i), meshes(idx))
       meshes(idx)%template_id = idx
       is_final = index(mesh_files(i), 'FIN/') > 0
       template_is_final(idx) = is_final
       template_loaded(idx) = .true.
    end do

    if (has_failure) then
       write(*, '(A)') 'Template catalogue contains unmatched entries; aborting.'
       stop 1
    end if

    write(*, '(A)', advance='no') 'Templates '
    do i = 1, mesh_file_count
       write(buf, '(I0)') i
       write(*, '(A)', advance='no') '[' // pad_field(buf) // ']'
    end do
    write(*, '(A)') ''
    write(*, '(A)', advance='no') 'Pattern   '
    do i = 1, mesh_file_count
       write(*, '(A)', advance='no') '[' // pattern_block(i) // ']'
    end do
    write(*, '(A)') ''
    write(*, '(A)', advance='no') 'Status    '
    do i = 1, mesh_file_count
       write(*, '(A)', advance='no') '[' // pad_field(merge('FIN','INT', template_is_final(i))) // ']'
    end do
    write(*, '(A)') ''
    write(*, '(A)', advance='no') 'Available '
    do i = 1, mesh_file_count
       write(buf, '(A)') merge('YES', 'NO ', template_loaded(i))
       write(*, '(A)', advance='no') '[' // pad_field(buf) // ']'
    end do
    write(*, '(A)') ''
    deallocate(template_loaded)
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

    pure function pattern_block(idx) result(block)
      integer, intent(in) :: idx
      character(len=8) :: block
      integer :: bit

      block = '--------'
      do bit = 1, min(8, size(templates, 1))
         if (templates(bit, idx)) block(bit:bit) = 'X'
      end do
    end function pattern_block
  end subroutine load_all_meshes

  subroutine release_all_meshes()
    integer :: i

    if (.not.allocated(meshes)) return
    do i = 1, size(meshes)
       call clear_mesh(meshes(i))
    end do
    deallocate(meshes)
  end subroutine release_all_meshes

  subroutine load_target_mesh(filename)
    character(len=*), intent(in) :: filename

    call clear_mesh(target_mesh)
    call read_mesh_file(filename, target_mesh)
    call update_cylindrical_radius(target_mesh)
    call load_setup_flags(filename)
  end subroutine load_target_mesh

  subroutine release_target_mesh()
    call clear_mesh(target_mesh)
  end subroutine release_target_mesh

  subroutine update_cylindrical_radius(mesh)
    type(mesh_type), intent(in) :: mesh
    integer :: ivt
    real(rk) :: radius, max_radius

    if (.not.allocated(mesh%coor)) then
       cylindrical_outer_radius = 0.0_rk
       return
    end if

    max_radius = 0.0_rk
    do ivt = 1, size(mesh%coor, 2)
       radius = sqrt(mesh%coor(1, ivt)**2 + mesh%coor(2, ivt)**2)
       if (radius > max_radius) max_radius = radius
    end do
    cylindrical_outer_radius = max_radius
  end subroutine update_cylindrical_radius

  subroutine read_mesh_file(filename, mesh)
    character(len=*), intent(in) :: filename
    type(mesh_type), intent(inout) :: mesh

    integer :: unit, ios, iel, ivt
    character(len=256) :: line
    integer :: nel, nvt, nbct, nve, nee, nae
    logical :: has_file

    inquire(file=trim(filename), exist=has_file)
    if (.not.has_file) then
       write(*, '(A)') 'Mesh file not found: ' // trim(filename)
       call clear_mesh(mesh)
       return
    end if

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*, '(A,I0)') 'Failed to open ' // trim(filename) // ' IOSTAT=', ios
      call clear_mesh(mesh)
      return
    end if

    read(unit, '(A)', iostat=ios) line
    read(unit, '(A)', iostat=ios) line
    read(unit, *, iostat=ios) nel, nvt, nbct, nve, nee, nae
    if (ios /= 0) then
      call abort_read('Failed to read header from ', filename, unit, mesh)
      return
    end if

    mesh%nel  = nel
    mesh%nvt  = nvt
    mesh%nbct = nbct
    mesh%nve  = nve
    mesh%nee  = nee
    mesh%nae  = nae

    read(unit, '(A)', iostat=ios) line  ! expect DCORVG tag
    if (ios /= 0) then
      call abort_read('Failed to read coordinate tag from ', filename, unit, mesh)
      return
    end if

    allocate(mesh%coor(3, nvt))
    do ivt = 1, nvt
       read(unit, *, iostat=ios) mesh%coor(1, ivt), mesh%coor(2, ivt), mesh%coor(3, ivt)
       if (ios /= 0) then
          call abort_read('Failed to read coordinates from ', filename, unit, mesh)
          return
       end if
    end do

    read(unit, '(A)', iostat=ios) line  ! expect KVERT tag
    if (ios /= 0) then
      call abort_read('Failed to read connectivity tag from ', filename, unit, mesh)
      return
    end if

    allocate(mesh%kvert(8, nel))
    do iel = 1, nel
       read(unit, *, iostat=ios) mesh%kvert(:, iel)
       if (ios /= 0) then
          call abort_read('Failed to read connectivities from ', filename, unit, mesh)
          return
       end if
    end do

    read(unit, '(A)', iostat=ios) line  ! expect KNPR tag
    if (ios /= 0) then
      call abort_read('Failed to read nodal property tag from ', filename, unit, mesh)
      return
    end if

    allocate(mesh%knpr(nvt))
    do ivt = 1, nvt
       read(unit, *, iostat=ios) mesh%knpr(ivt)
       if (ios /= 0) then
          call abort_read('Failed to read nodal properties from ', filename, unit, mesh)
          return
       end if
    end do

    close(unit)
  end subroutine read_mesh_file

  subroutine abort_read(message, filename, unit, mesh)
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: unit
    type(mesh_type), intent(inout) :: mesh

    write(*, '(A)') trim(message) // trim(filename)
    close(unit)
    call clear_mesh(mesh)
  end subroutine abort_read

  subroutine build_clean_output_path(base_path, clean_path)
    character(len=*), intent(in) :: base_path
    character(len=*), intent(out) :: clean_path
    character(len=1024) :: work, prefix, suffix
    integer :: len_in, i, last_sep, last_dot

    clean_path = ''
    work = trim(base_path)
    len_in = len_trim(work)
    if (len_in <= 0) then
       clean_path = 'refined_mesh_clean.vtu'
       return
    end if

    last_sep = 0
    last_dot = 0
    do i = 1, len_in
       select case (work(i:i))
       case ('/', '\')
          last_sep = i
       case ('.')
          last_dot = i
       end select
    end do
    if (last_dot <= last_sep) last_dot = 0

    prefix = ''
    suffix = ''
    if (last_dot > 0) then
       if (last_dot > 1) prefix = work(1:last_dot-1)
       suffix = work(last_dot:len_in)
       clean_path = trim(prefix) // '_clean' // trim(suffix)
    else
       clean_path = work(1:len_in) // '_clean'
    end if
  end subroutine build_clean_output_path

  subroutine build_refined_clean_output_path(base_path, clean_path)
    character(len=*), intent(in) :: base_path
    character(len=*), intent(out) :: clean_path
    character(len=1024) :: work, prefix, suffix
    integer :: len_in, i, last_sep, last_dot

    clean_path = ''
    work = trim(base_path)
    len_in = len_trim(work)
    if (len_in <= 0) then
       clean_path = 'refined_mesh_refined_clean.vtu'
       return
    end if

    last_sep = 0
    last_dot = 0
    do i = 1, len_in
       select case (work(i:i))
       case ('/', '\')
          last_sep = i
       case ('.')
          last_dot = i
       end select
    end do
    if (last_dot <= last_sep) last_dot = 0

    prefix = ''
    suffix = ''
    if (last_dot > 0) then
      if (last_dot > 1) prefix = work(1:last_dot-1)
      suffix = work(last_dot:len_in)
      clean_path = trim(prefix) // '_refined_clean' // trim(suffix)
    else
      clean_path = work(1:len_in) // '_refined_clean'
    end if
  end subroutine build_refined_clean_output_path

  subroutine build_level_output_path(base_path, level, level_path)
    character(len=*), intent(in) :: base_path
    integer, intent(in) :: level
    character(len=*), intent(out) :: level_path
    character(len=1024) :: work, prefix, suffix
    character(len=32) :: level_str
    integer :: len_in, i, last_sep, last_dot

    level_path = ''
    work = trim(base_path)
    len_in = len_trim(work)
    write(level_str, '(I0)') level
    if (len_in <= 0) then
       level_path = 'refined_mesh_lvl' // trim(level_str) // '.vtu'
       return
    end if

    last_sep = 0
    last_dot = 0
    do i = 1, len_in
       select case (work(i:i))
       case ('/', '\')
          last_sep = i
       case ('.')
          last_dot = i
       end select
    end do
    if (last_dot <= last_sep) last_dot = 0

    prefix = ''
    suffix = ''
    if (last_dot > 0) then
       if (last_dot > 1) prefix = work(1:last_dot-1)
       suffix = work(last_dot:len_in)
       level_path = trim(prefix) // '_lvl' // trim(level_str) // trim(suffix)
    else
       level_path = work(1:len_in) // '_lvl' // trim(level_str)
    end if
  end subroutine build_level_output_path

  subroutine write_patch_group_vtu(patch_groups, filename)
    type(element_patch_group), intent(in) :: patch_groups(:)
    character(len=*), intent(in) :: filename
    integer :: total_points, total_cells
    integer :: i, j, k, point_offset, cell_count, conn_idx, elem_idx
    real(rk), allocatable :: points(:,:)
    integer, allocatable   :: point_knpr(:)
    integer, allocatable   :: cell_conn(:)
    integer, allocatable   :: offsets(:)
    integer, allocatable   :: cell_types(:)
    integer, allocatable   :: cell_monitor(:)
    integer :: unit, ios, knpr_value
    character(len=32) :: num_points_str, num_cells_str

    if (size(patch_groups) <= 0) then
       write(*, '(A)') 'No patch groups available; VTU file not written.'
       return
    end if

    total_points = 0
    total_cells  = 0
    do elem_idx = 1, size(patch_groups)
       if (.not.allocated(patch_groups(elem_idx)%patchlist)) cycle
       do i = 1, patch_groups(elem_idx)%count
          associate(patch => patch_groups(elem_idx)%patchlist(i))
             if (patch%n_vert > 0 .and. patch%n_elem > 0) then
                total_points = total_points + patch%n_vert
                total_cells  = total_cells  + patch%n_elem
             end if
          end associate
       end do
    end do

    if (total_points == 0 .or. total_cells == 0) then
      write(*, '(A)') 'No populated element patches; VTU file not written.'
      return
    end if

    allocate(points(3, total_points))
    allocate(point_knpr(total_points))
    allocate(cell_conn(8 * total_cells))
    allocate(offsets(total_cells))
    allocate(cell_types(total_cells))
    allocate(cell_monitor(total_cells))

    point_knpr = 0
    cell_conn  = 0
    offsets    = 0
    cell_types = 12
    cell_monitor = 0

    point_offset = 0
    cell_count   = 0
    conn_idx     = 0
    do elem_idx = 1, size(patch_groups)
       if (.not.allocated(patch_groups(elem_idx)%patchlist)) cycle
       do i = 1, patch_groups(elem_idx)%count
          associate(patch => patch_groups(elem_idx)%patchlist(i))
             if (patch%n_vert <= 0 .or. patch%n_elem <= 0) cycle
             if (.not.allocated(patch%global_coor)) cycle
             if (.not.allocated(patch%kvert)) cycle

             do j = 1, patch%n_vert
                points(:, point_offset + j) = patch%global_coor(:, j)
                if (allocated(patch%knpr) .and. j <= size(patch%knpr)) then
                   knpr_value = patch%knpr(j)
                else
                   knpr_value = 0
                end if
                point_knpr(point_offset + j) = knpr_value
             end do

             do j = 1, patch%n_elem
                cell_count = cell_count + 1
                if (allocated(patch%monitor)) then
                   if (j <= size(patch%monitor)) then
                      cell_monitor(cell_count) = patch%monitor(j)
                   else
                      cell_monitor(cell_count) = 0
                   end if
                else
                   cell_monitor(cell_count) = 0
                end if
                offsets(cell_count) = cell_count * 8
                do k = 1, 8
                   conn_idx = conn_idx + 1
                   cell_conn(conn_idx) = point_offset + patch%kvert(k, j) - 1
                end do
             end do

             point_offset = point_offset + patch%n_vert
          end associate
       end do
    end do

    total_cells = cell_count
    total_points = point_offset

    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
       write(*, '(A,I0)') 'Failed to open VTU file, IOSTAT=', ios
       deallocate(points, point_knpr, cell_conn, offsets, cell_types)
       return
    end if

    write(num_points_str, '(I0)') total_points
    write(num_cells_str, '(I0)') total_cells

    write(unit, '(A)') '<?xml version="1.0"?>'
    write(unit, '(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(unit, '(A)') '  <UnstructuredGrid>'
    write(unit, '(A)') '    <Piece NumberOfPoints="' // trim(num_points_str) // &
         '" NumberOfCells="' // trim(num_cells_str) // '">'

    write(unit, '(A)') '      <PointData Scalars="KNPR">'
    write(unit, '(A)') '        <DataArray type="Int32" Name="KNPR" format="ascii">'
    do i = 1, total_points, 6
       write(unit, '(6(I12))') point_knpr(i:min(i+5, total_points))
    end do
    write(unit, '(A)') '        </DataArray>'
    write(unit, '(A)') '      </PointData>'
    write(unit, '(A)') '      <CellData>'
    write(unit, '(A)') '        <DataArray type="Int32" Name="monitor" format="ascii">'
    do i = 1, total_cells, 8
       write(unit, '(8(I12))') cell_monitor(i:min(i+7, total_cells))
    end do
    write(unit, '(A)') '        </DataArray>'
    write(unit, '(A)') '      </CellData>'

    write(unit, '(A)') '      <Points>'
    write(unit, '(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    do i = 1, total_points
       write(unit, '(3(ES24.16,1X))') points(1, i), points(2, i), points(3, i)
    end do
    write(unit, '(A)') '        </DataArray>'
    write(unit, '(A)') '      </Points>'

    write(unit, '(A)') '      <Cells>'
    write(unit, '(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
    do i = 1, total_cells
       write(unit, '(8(I12))') cell_conn(8*(i-1)+1:8*i)
    end do
    write(unit, '(A)') '        </DataArray>'

    write(unit, '(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
    do i = 1, total_cells
       write(unit, '(I12)') offsets(i)
    end do
    write(unit, '(A)') '        </DataArray>'

    write(unit, '(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
    do i = 1, total_cells
       write(unit, '(I12)') cell_types(i)
    end do
    write(unit, '(A)') '        </DataArray>'
    write(unit, '(A)') '      </Cells>'

    write(unit, '(A)') '    </Piece>'
    write(unit, '(A)') '  </UnstructuredGrid>'
    write(unit, '(A)') '</VTKFile>'
    close(unit)

    deallocate(points, point_knpr, cell_conn, offsets, cell_types, cell_monitor)
    write(*, '(A)') 'VTU file written: ' // trim(filename)
  end subroutine write_patch_group_vtu

  subroutine write_refined_clean_tri(mesh, reference_mesh, working_folder)
    type(mesh_type), intent(in) :: mesh
    type(mesh_type), intent(in) :: reference_mesh
    character(len=*), intent(in) :: working_folder

    integer :: unit, ios, iel, ivt
    integer :: nbct, nve, nee, nae
    character(len=1024) :: folder_path, file_path
    real(rk) :: scaled_point(3)

    if (mesh%nel <= 0 .or. mesh%nvt <= 0) then
       write(*, '(A)') 'Refined mesh is empty; TRI export skipped.'
       return
    end if
    if (.not.allocated(mesh%coor) .or. .not.allocated(mesh%kvert)) then
       write(*, '(A)') 'Refined mesh missing coordinates or connectivity; TRI export skipped.'
       return
    end if

    folder_path = trim(adjustl(working_folder))
    if (len_trim(folder_path) == 0) folder_path = '.'
    folder_path = trim(folder_path) // '/meshDir_BU'
    if (.not.ensure_directory(folder_path)) then
       write(*, '(A)') 'Failed to create output folder for TRI export: ' // trim(folder_path)
       return
    end if

    file_path = trim(folder_path) // '/Merged_Mesh.tri'

    nbct = reference_mesh%nbct
    nve  = reference_mesh%nve
    nee  = reference_mesh%nee
    nae  = reference_mesh%nae

    open(newunit=unit, file=trim(file_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
       write(*, '(A,I0)') 'Failed to open TRI file for write, IOSTAT=', ios
       return
    end if

    write(unit, '(A)') 'Coarse mesh exported by DeViSoR TRI3D exporter'
    write(unit, '(A)') 'Parametrisierung PARXC, PARYC, TMAXC'
    write(unit, '(6(I12,1X),5X,A)') mesh%nel, mesh%nvt, nbct, nve, nee, nae, '    NEL,NVT,NBCT,NVE,NEE,NAE'
    write(unit, '(A)') 'DCORVG'
    do ivt = 1, mesh%nvt
       scaled_point = scaling_factor_TRI_output * mesh%coor(:, ivt)
       write(unit, '(3(1X,ES24.16))') scaled_point(1), scaled_point(2), scaled_point(3)
    end do
    write(unit, '(A)') 'KVERT'
    do iel = 1, mesh%nel
       write(unit, '(8(I12))') mesh%kvert(:, iel)
    end do
    write(unit, '(A)') 'KNPR'
    do ivt = 1, mesh%nvt
       write(unit, '(I12)') 0
    end do

    close(unit)
    write(*, '(A)') 'TRI file written: ' // trim(file_path)
  end subroutine write_refined_clean_tri

  subroutine load_setup_flags(mesh_filename)
    character(len=*), intent(in) :: mesh_filename
    character(len=1024) :: setup_path, line, key, mesh_dir, working_dir
    integer :: unit, ios
    logical :: in_section, has_file

    mesh_dir = trim(extract_parent_folder(mesh_filename))
    working_dir = trim(extract_parent_folder(mesh_dir))
    setup_path = trim(working_dir) // '/setup.e3d'
    inquire(file=trim(setup_path), exist=has_file)
    if (.not.has_file) then
       write(*, '(A)') 'setup.e3d not found; defaulting to box transformation.'
       apply_cylindric_transform = .false.
       return
    end if

    open(newunit=unit, file=trim(setup_path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*, '(A)') 'Unable to open setup.e3d; defaulting to box transformation.'
       apply_cylindric_transform = .false.
       return
    end if
    in_section = .false.
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (line(1:1) == '[') then
          if (index(line, 'E3DGeometryData/Preprocessing') > 0) then
             in_section = .true.
          else
             in_section = .false.
          end if
          cycle
       end if
       if (.not.in_section) cycle
       key = adjustl(line)
       if (index(key, 'HexMesher=') > 0) then
          if (index(key, 'HollowCylinder') > 0 .or. index(key, 'FullCylinder') > 0) then
             apply_cylindric_transform = .true.
             write(*, '(A)') 'Cylindrical transformation enabled based on setup.e3d (HexMesher requests cylindrical mesh).'
          else
             apply_cylindric_transform = .false.
             write(*, '(A)') 'Box transformation selected (HexMesher setting is not cylindrical).'
          end if
          exit
       end if
    end do
    close(unit)
  end subroutine load_setup_flags

  logical function ensure_directory(path)
    character(len=*), intent(in) :: path
    character(len=2048) :: command
    integer :: exit_status

    command = 'mkdir -p "' // trim(path) // '"'
    call execute_command_line(trim(command), exitstat=exit_status)
    ensure_directory = (exit_status == 0)
  end function ensure_directory

  character(len=1024) pure function extract_parent_folder(path)
    character(len=*), intent(in) :: path
    integer :: i, last_sep

    extract_parent_folder = ''
    last_sep = 0
    do i = 1, len_trim(path)
       if (path(i:i) == '/' .or. path(i:i) == '\') last_sep = i
    end do
    if (last_sep <= 1) then
       extract_parent_folder = '.'
    else
       extract_parent_folder = path(1:last_sep-1)
    end if
  end function extract_parent_folder

end module inout_mod
