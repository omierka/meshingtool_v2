module hex_io
  use tri_tet_intersection, only: dp
  implicit none
  private
  public :: read_hex_mesh
  public :: read_tri_mesh
  public :: write_hex_mesh_vtu
  public :: write_hex_intersection_vtu
  public :: write_pvtu_reference

contains

  subroutine read_hex_mesh(filename, dcorvg, kvert, nvt, nel)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: dcorvg(:, :)
    integer, allocatable, intent(out) :: kvert(:, :)
    integer, intent(out) :: nvt, nel

    integer :: unit, ios
    character(len=256) :: line
    integer :: nbct, nve, nee, nae
    integer :: i

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Failed to open mesh file ', trim(filename)
      stop 1
    end if

    read(unit,'(A)', iostat=ios) line
    read(unit,'(A)', iostat=ios) line
    read(unit,'(A)', iostat=ios) line
    if (ios /= 0) then
      write(*,*) 'Mesh header too short in ', trim(filename)
      stop 1
    end if
    read(line,*) nel, nvt, nbct, nve, nee, nae

    read(unit,'(A)', iostat=ios) line
    if (ios /= 0) then
      write(*,*) 'Missing DCORVG label in ', trim(filename)
      stop 1
    end if

    allocate(dcorvg(3, nvt))
    do i = 1, nvt
      read(unit,*, iostat=ios) dcorvg(1,i), dcorvg(2,i), dcorvg(3,i)
      if (ios /= 0) then
        write(*,*) 'Failed reading coordinates at index ', i
        stop 1
      end if
    end do

    do
      read(unit,'(A)', iostat=ios) line
      if (ios /= 0) then
        write(*,*) 'Missing KVERT section in ', trim(filename)
        stop 1
      end if
      if (len_trim(line) == 0) cycle
      line = adjustl(line)
      if (len_trim(line) >= 5) then
        if (line(1:5) == 'KVERT') exit
      end if
    end do

    allocate(kvert(8, nel))
    do i = 1, nel
      read(unit,*, iostat=ios) kvert(:, i)
      if (ios /= 0) then
        write(*,*) 'Failed reading connectivity at element ', i
        stop 1
      end if
    end do

    close(unit)
  end subroutine read_hex_mesh

  subroutine read_tri_mesh(filename, dcorvg, kvert, nvt, ntri, MIS_diameter)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: dcorvg(:, :)
    integer, allocatable, intent(out) :: kvert(:, :)
    integer, intent(out) :: nvt, ntri
    real(dp), allocatable, intent(out) :: MIS_diameter(:)

    character(len=256) :: fname_trim
    character(len=4) :: ext
    integer :: lenf

    fname_trim = trim(filename)
    lenf = len_trim(fname_trim)
    ext = ' '
    if (lenf >= 4) then
      ext = fname_trim(lenf-3:lenf)
      call to_lowercase(ext)
    end if

    if (ext == '.vtu') then
      call read_tri_mesh_vtu(fname_trim, dcorvg, kvert, nvt, ntri, MIS_diameter)
    else
      call read_tri_mesh_off(fname_trim, dcorvg, kvert, nvt, ntri, MIS_diameter)
    end if

  contains

    subroutine to_lowercase(str)
      character(len=*), intent(inout) :: str
      integer :: i, ia
      do i = 1, len(str)
        ia = iachar(str(i:i))
        if (ia >= iachar('A') .and. ia <= iachar('Z')) then
          str(i:i) = achar(ia + (iachar('a') - iachar('A')))
        end if
      end do
    end subroutine to_lowercase

  end subroutine read_tri_mesh

  subroutine read_tri_mesh_off(filename, dcorvg, kvert, nvt, ntri, MIS_diameter)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: dcorvg(:, :)
    integer, allocatable, intent(out) :: kvert(:, :)
    integer, intent(out) :: nvt, ntri
    real(dp), allocatable, intent(out) :: MIS_diameter(:)

    integer :: unit, ios
    character(len=256) :: line
    integer :: i, nignore
    integer :: nv_face

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Failed to open triangle mesh file ', trim(filename)
      stop 1
    end if

    read(unit,'(A)', iostat=ios) line
    if (ios /= 0) then
      write(*,*) 'Triangle mesh header missing in ', trim(filename)
      stop 1
    end if

    read(unit,'(A)', iostat=ios) line
    if (ios /= 0) then
      write(*,*) 'Triangle mesh count line missing in ', trim(filename)
      stop 1
    end if
    read(line,*) nvt, ntri, nignore
    if (nvt <= 0 .or. ntri <= 0) then
      write(*,*) 'Invalid counts in triangle mesh header'
      stop 1
    end if

    allocate(dcorvg(3, nvt))
    do i = 1, nvt
      read(unit,*, iostat=ios) dcorvg(1,i), dcorvg(2,i), dcorvg(3,i)
      if (ios /= 0) then
        write(*,*) 'Failed reading triangle coordinate at index ', i
        stop 1
      end if
    end do

    allocate(kvert(3, ntri))
    do i = 1, ntri
      read(unit,*, iostat=ios) nv_face, kvert(1, i), kvert(2, i), kvert(3, i)
      if (ios /= 0) then
        write(*,*) 'Failed reading triangle connectivity at element ', i
        stop 1
      end if
      if (nv_face /= 3) then
        write(*,*) 'Encountered non-triangular face in ', trim(filename)
        stop 1
      end if
      kvert(:, i) = kvert(:, i) + 1
    end do
    close(unit)

    allocate(MIS_diameter(ntri))
    MIS_diameter = 0.0_dp
  end subroutine read_tri_mesh_off

  subroutine read_tri_mesh_vtu(filename, dcorvg, kvert, nvt, ntri, MIS_diameter)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: dcorvg(:, :)
    integer, allocatable, intent(out) :: kvert(:, :)
    integer, intent(out) :: nvt, ntri
    real(dp), allocatable, intent(out) :: MIS_diameter(:)

    integer :: unit, ios
    character(len=2048) :: line
    character(len=2048) :: trimmed
    integer :: nCells
    logical :: have_piece
    integer, parameter :: section_none = 0, section_points = 1, section_conn = 2, &
         section_offsets = 3, section_types = 4, section_mis = 5
    integer :: active_section
    real(dp), allocatable :: point_values(:)
    integer, allocatable :: conn_values(:)
    integer, allocatable :: offsets(:)
    integer, allocatable :: cell_types(:)
    real(dp), allocatable :: cell_mis(:)
    integer :: point_count, conn_count, offset_count, type_count, mis_count
    integer :: i, prev_offset, tri_idx, cell_vertices, start_idx
    integer :: last_offset

    nvt = -1
    nCells = -1
    have_piece = .false.
    active_section = section_none
    point_count = 0
    conn_count = 0
    offset_count = 0
    type_count = 0
    mis_count = 0

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Failed to open triangle VTU file ', trim(filename)
      stop 1
    end if

    do
      read(unit,'(A)', iostat=ios) line
      if (ios /= 0) exit
      trimmed = adjustl(line)
      if (.not. have_piece) then
        if (index(trimmed, '<Piece') > 0) then
          call extract_piece_counts(trimmed, nvt, nCells)
          if (nvt <= 0 .or. nCells <= 0) then
            write(*,*) 'Invalid Piece header in ', trim(filename)
            stop 1
          end if
          allocate(point_values(3*nvt))
          allocate(conn_values(max(3*nCells, 1)))
          allocate(offsets(nCells))
          allocate(cell_types(nCells))
          allocate(cell_mis(nCells))
          have_piece = .true.
        end if
        cycle
      end if

      if (index(trimmed, '<DataArray') > 0) then
        active_section = identify_section(trimmed)
        call process_data_line(trimmed, active_section)
        if (index(trimmed, '</DataArray>') > 0) active_section = section_none
        cycle
      end if

      if (active_section /= section_none) then
        call process_data_line(trimmed, active_section)
        if (index(trimmed, '</DataArray>') > 0) active_section = section_none
      end if
    end do

    close(unit)

    if (.not. have_piece) then
      write(*,*) 'No <Piece> section found in ', trim(filename)
      stop 1
    end if
    if (point_count /= 3*nvt) then
      write(*,*) 'Unexpected number of point coordinates in ', trim(filename)
      stop 1
    end if
    !if (offset_count /= nCells .or. type_count /= nCells .or. mis_count /= nCells) then
    !  write(*,*) 'Mismatch in cell data counts in ', trim(filename),offset_count , type_count  , mis_count, nCells
    !  stop 1
    !end if

    last_offset = offsets(nCells)
    if (conn_count < last_offset) then
      write(*,*) 'Connectivity array shorter than offsets in ', trim(filename)
      stop 1
    end if

    ntri = count(cell_types == 5)
    if (ntri <= 0) then
      write(*,*) 'No triangle cells found in ', trim(filename)
      stop 1
    end if

    allocate(dcorvg(3, nvt))
    do i = 1, nvt
      dcorvg(1,i) = point_values(3*(i-1)+1)
      dcorvg(2,i) = point_values(3*(i-1)+2)
      dcorvg(3,i) = point_values(3*(i-1)+3)
    end do

    allocate(kvert(3, ntri))
    allocate(MIS_diameter(ntri))
    prev_offset = 0
    tri_idx = 0
    do i = 1, nCells
      cell_vertices = offsets(i) - prev_offset
      start_idx = prev_offset + 1
      if (cell_types(i) == 5) then
        if (cell_vertices /= 3) then
          write(*,*) 'Triangle cell with ', cell_vertices, ' nodes encountered'
          stop 1
        end if
        tri_idx = tri_idx + 1
        kvert(:, tri_idx) = conn_values(start_idx:start_idx+2) + 1
        MIS_diameter(tri_idx) = cell_mis(i)
      end if
      prev_offset = offsets(i)
    end do

    deallocate(point_values, conn_values, offsets, cell_types, cell_mis)

  contains

    subroutine extract_piece_counts(text, points, cells)
      character(len=*), intent(in) :: text
      integer, intent(out) :: points, cells
      points = read_attribute_int(text, 'NumberOfPoints')
      cells = read_attribute_int(text, 'NumberOfCells')
    end subroutine extract_piece_counts

    integer function read_attribute_int(text, attr)
      character(len=*), intent(in) :: text
      character(len=*), intent(in) :: attr
      integer :: pos_attr, start_pos, end_pos, ios_local, len_text, copy_len
      character(len=64) :: buffer
      read_attribute_int = -1
      len_text = len_trim(text)
      pos_attr = index(text, trim(attr)//'="')
      if (pos_attr <= 0) return
      start_pos = pos_attr + len_trim(attr) + 2
      end_pos = start_pos
      do while (end_pos <= len_text)
        if (text(end_pos:end_pos) == '"') exit
        end_pos = end_pos + 1
      end do
      if (end_pos > len_text) return
      if (end_pos <= start_pos) return
      copy_len = min(end_pos - start_pos, len(buffer))
      buffer = ' '
      buffer(1:copy_len) = text(start_pos:start_pos+copy_len-1)
      read(buffer(1:copy_len),*,iostat=ios_local) read_attribute_int
      if (ios_local /= 0) read_attribute_int = -1
    end function read_attribute_int

    integer function identify_section(text)
      character(len=*), intent(in) :: text
      if (index(text, 'Name="Points"') > 0) then
        identify_section = section_points
      else if (index(text, 'Name="connectivity"') > 0) then
        identify_section = section_conn
      else if (index(text, 'Name="offsets"') > 0) then
        identify_section = section_offsets
      else if (index(text, 'Name="types"') > 0) then
        identify_section = section_types
      else if (index(text, 'Name="Monitor"') > 0) then
        identify_section = section_mis
      else
        identify_section = section_none
      end if
    end function identify_section

    subroutine process_data_line(text, section_id)
      character(len=*), intent(in) :: text
      integer, intent(in) :: section_id
      integer :: gtpos, endpos
      character(len=:), allocatable :: payload

      if (section_id == section_none) return
      gtpos = index(text, '>')
      if (gtpos > 0) then
        payload = text(gtpos+1:)
      else
        payload = text
      end if

      endpos = index(payload, '</DataArray>')
      if (endpos > 0) then
        payload = payload(:endpos-1)
      end if

      if (len_trim(payload) == 0) then
        if (index(text, '</DataArray>') > 0) active_section = section_none
        return
      end if

      select case (section_id)
      case (section_points)
        call append_reals(payload, point_values, point_count, 3*nvt)
      case (section_conn)
        call append_ints(payload, conn_values, conn_count, -1, .true.)
      case (section_offsets)
        call append_ints(payload, offsets, offset_count, nCells)
      case (section_types)
        call append_ints(payload, cell_types, type_count, nCells)
      case (section_mis)
        call append_reals(payload, cell_mis, mis_count, nCells)
      end select

      if (index(text, '</DataArray>') > 0) active_section = section_none
    end subroutine process_data_line

    subroutine append_reals(text, buffer, count, expected)
      character(len=*), intent(in) :: text
      real(dp), intent(inout) :: buffer(:)
      integer, intent(inout) :: count
      integer, intent(in) :: expected
      integer :: len_text, i, token_len
      character(len=64) :: token
      character :: ch

      len_text = len_trim(text)
      if (len_text <= 0) return
      token_len = 0
      token = ' '
      do i = 1, len_text
        ch = text(i:i)
        if (is_real_char(ch)) then
          token_len = token_len + 1
          if (token_len <= len(token)) token(token_len:token_len) = ch
        else
          if (token_len > 0) then
            call store_real_value(token, token_len, buffer, count, expected)
            token_len = 0
            token = ' '
          end if
        end if
      end do
      if (token_len > 0) call store_real_value(token, token_len, buffer, count, expected)
    end subroutine append_reals

    subroutine append_ints(text, buffer, count, expected, allow_resize)
      character(len=*), intent(in) :: text
      integer, allocatable, intent(inout) :: buffer(:)
      integer, intent(inout) :: count
      integer, intent(in) :: expected
      logical, intent(in), optional :: allow_resize
      integer :: len_text, i, token_len
      character(len=64) :: token
      character :: ch
      logical :: do_resize

      len_text = len_trim(text)
      if (len_text <= 0) return
      token_len = 0
      token = ' '
      do_resize = .false.
      if (present(allow_resize)) do_resize = allow_resize
      do i = 1, len_text
        ch = text(i:i)
        if (is_int_char(ch)) then
          token_len = token_len + 1
          if (token_len <= len(token)) token(token_len:token_len) = ch
        else
          if (token_len > 0) then
            call store_int_value(token, token_len, buffer, count, expected, do_resize)
            token_len = 0
            token = ' '
          end if
        end if
      end do
      if (token_len > 0) call store_int_value(token, token_len, buffer, count, expected, do_resize)
    end subroutine append_ints

    subroutine store_real_value(token, token_len, buffer, count, expected)
      character(len=*), intent(in) :: token
      integer, intent(in) :: token_len
      real(dp), intent(inout) :: buffer(:)
      integer, intent(inout) :: count
      integer, intent(in) :: expected
      integer :: ios_local
      real(dp) :: value

      if (token_len <= 0) return
      read(token(1:token_len),*,iostat=ios_local) value
      if (ios_local /= 0) return
      if (count >= size(buffer) .or. count >= expected) then
        write(*,*) 'Exceeded real value allocation while reading ', trim(filename)
        stop 1
      end if
      count = count + 1
      buffer(count) = value
    end subroutine store_real_value

    subroutine store_int_value(token, token_len, buffer, count, expected, allow_resize)
      character(len=*), intent(in) :: token
      integer, intent(in) :: token_len
      integer, allocatable, intent(inout) :: buffer(:)
      integer, intent(inout) :: count
      integer, intent(in) :: expected
      logical, intent(in) :: allow_resize
      integer :: ios_local
      integer :: value

      if (token_len <= 0) return
      read(token(1:token_len),*,iostat=ios_local) value
      if (ios_local /= 0) return
      call ensure_int_capacity(buffer, count, expected, allow_resize)
      count = count + 1
      buffer(count) = value
    end subroutine store_int_value

    subroutine ensure_int_capacity(buffer, count, expected, allow_resize)
      integer, allocatable, intent(inout) :: buffer(:)
      integer, intent(in) :: count
      integer, intent(in) :: expected
      logical, intent(in) :: allow_resize
      integer :: new_size, old_size
      integer, allocatable :: tmp(:)

      if (expected >= 0 .and. count >= expected) then
        write(*,*) 'Exceeded integer value allocation while reading ', trim(filename)
        stop 1
      end if
      if (count < size(buffer)) return
      if (.not. allow_resize) then
        write(*,*) 'Exceeded integer value allocation while reading ', trim(filename)
        stop 1
      end if
      old_size = size(buffer)
      new_size = max(count + 1, max(2*old_size, 1))
      allocate(tmp(new_size))
      tmp(1:old_size) = buffer
      deallocate(buffer)
      call move_alloc(tmp, buffer)
    end subroutine ensure_int_capacity

    logical function is_real_char(ch)
      character, intent(in) :: ch
      is_real_char = (ch >= '0' .and. ch <= '9') .or. ch == '-' .or. ch == '+' .or. &
           ch == '.' .or. ch == 'e' .or. ch == 'E' .or. ch == 'd' .or. ch == 'D'
    end function is_real_char

    logical function is_int_char(ch)
      character, intent(in) :: ch
      is_int_char = (ch >= '0' .and. ch <= '9') .or. ch == '-' .or. ch == '+'
    end function is_int_char

  end subroutine read_tri_mesh_vtu

  subroutine write_hex_intersection_vtu(filename_base, tet_points, tet_conn, &
       int_points, int_conn, int_offsets, int_types, int_hex_ids, tet_file_out, int_file_out)
    character(len=*), intent(in) :: filename_base
    real(dp), intent(in) :: tet_points(:, :)
    integer, intent(in) :: tet_conn(:, :)
    real(dp), intent(in) :: int_points(:, :)
    integer, intent(in) :: int_conn(:)
    integer, intent(in) :: int_offsets(:)
    integer, intent(in) :: int_types(:)
    integer, intent(in) :: int_hex_ids(:)
    character(len=*), intent(out), optional :: tet_file_out, int_file_out

    integer :: total_tets, total_int_points, total_int_cells, total_int_conn
    integer, allocatable :: tet_conn_flat(:), tet_offsets(:), tet_types(:)
    character(len=512) :: filename_trim, dir_path, base_name
    character(len=512) :: tet_file_name, int_file_name
    character(len=512) :: tet_file_path, int_file_path
    integer :: lenf, slashpos, ic, i

    total_tets = size(tet_conn, 2)
    total_int_points = size(int_points, 2)
    total_int_conn = size(int_conn)
    total_int_cells = size(int_offsets)

    allocate(tet_conn_flat(4*total_tets), tet_offsets(total_tets), tet_types(total_tets))
    tet_offsets = [(4*i, i=1,total_tets)]
    tet_types = 10
    if (total_tets > 0) then
      tet_conn_flat = reshape(tet_conn, [4*total_tets]) - 1
    end if

    filename_trim = trim(filename_base)
    dir_path = ''
    base_name = filename_trim
    lenf = len_trim(filename_trim)
    slashpos = 0
    do ic = lenf, 1, -1
      if (filename_trim(ic:ic) == '/' .or. filename_trim(ic:ic) == '\') then
        slashpos = ic
        exit
      end if
    end do
    if (slashpos > 0) then
      dir_path = filename_trim(:slashpos)
      base_name = filename_trim(slashpos+1:)
    end if

    tet_file_name = trim(base_name)//'_tets.vtu'
    int_file_name = trim(base_name)//'_intersections.vtu'
    tet_file_path = trim(dir_path)//trim(tet_file_name)
    int_file_path = trim(dir_path)//trim(int_file_name)

    call write_vtu_file(tet_file_path, tet_points, size(tet_points,2), tet_conn_flat, 4*total_tets, &
         tet_offsets, total_tets, tet_types, total_tets)
    call write_vtu_file(int_file_path, int_points, total_int_points, int_conn, total_int_conn, &
         int_offsets, total_int_cells, int_types, total_int_cells, int_hex_ids, "hex_id")

    if (present(tet_file_out)) tet_file_out = tet_file_path
    if (present(int_file_out)) int_file_out = int_file_path
    deallocate(tet_conn_flat, tet_offsets, tet_types)

  contains

    subroutine write_vtu_file(path, points, nPoints, conn, nConn, offsets, nOffsets, types, nTypes, cellField, fieldName)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: points(:, :)
      integer, intent(in) :: nPoints
      integer, intent(in) :: conn(:), nConn
      integer, intent(in) :: offsets(:), nOffsets
      integer, intent(in) :: types(:), nTypes
      integer, intent(in), optional :: cellField(:)
      character(len=*), intent(in), optional :: fieldName
      integer :: unitFile, iosFile, i

      open(newunit=unitFile, file=path, status='replace', action='write', iostat=iosFile)
      if (iosFile /= 0) then
        write(*,*) 'Failed to open file ', trim(path)
        return
      end if

      write(unitFile,'(A)') '<?xml version="1.0"?>'
      write(unitFile,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
      write(unitFile,'(A)') '  <UnstructuredGrid>'
      write(unitFile,'(A," NumberOfPoints=""",I0,""" NumberOfCells=""",I0,""">")') '    <Piece', nPoints, nOffsets
      write(unitFile,'(A)') '      <PointData/>'
      if (present(cellField)) then
        write(unitFile,'(A," Scalars=""",A,""">")') '      <CellData', trim(fieldName)
        write(unitFile,'(A," type=""Int32"" Name=""",A,""" format=""ascii"">")') '        <DataArray', trim(fieldName)
        write(unitFile,'(A,*(I0,1X))') '          ', cellField
        write(unitFile,'(A)') '        </DataArray>'
        write(unitFile,'(A)') '      </CellData>'
      else
        write(unitFile,'(A)') '      <CellData/>'
      end if
      write(unitFile,'(A)') '      <Points>'
      write(unitFile,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
      if (nPoints > 0) then
        do i = 1, nPoints
          write(unitFile,'(A,3(1X,ES23.15))') '          ', points(1, i), points(2, i), points(3, i)
        end do
      end if
      write(unitFile,'(A)') '        </DataArray>'
      write(unitFile,'(A)') '      </Points>'
      write(unitFile,'(A)') '      <Cells>'
      write(unitFile,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
      if (nConn > 0) then
        write(unitFile,'(A,*(I0,1X))') '          ', conn
      end if
      write(unitFile,'(A)') '        </DataArray>'
      write(unitFile,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
      if (nOffsets > 0) then
        write(unitFile,'(A,*(I0,1X))') '          ', offsets
      end if
      write(unitFile,'(A)') '        </DataArray>'
      write(unitFile,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
      if (nTypes > 0) then
        write(unitFile,'(A,*(I0,1X))') '          ', types
      end if
      write(unitFile,'(A)') '        </DataArray>'
      write(unitFile,'(A)') '      </Cells>'
      write(unitFile,'(A)') '    </Piece>'
      write(unitFile,'(A)') '  </UnstructuredGrid>'
      write(unitFile,'(A)') '</VTKFile>'
      close(unitFile)
    end subroutine write_vtu_file

  end subroutine write_hex_intersection_vtu

  subroutine write_hex_mesh_vtu(filename, points, connectivity, elem_indices, cell_field, field_name, &
       cell_field2, field_name2)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: points(:, :)
    integer, intent(in) :: connectivity(:, :)
    integer, intent(in), optional :: elem_indices(:)
    real(dp), intent(in), optional :: cell_field(:)
    character(len=*), intent(in), optional :: field_name
    real(dp), intent(in), optional :: cell_field2(:)
    character(len=*), intent(in), optional :: field_name2

    integer :: nel, nvt
    integer :: unit, ios
    integer, allocatable :: conn_flat(:), offsets(:), cellTypes(:)
    integer :: i, idx
    logical :: has_cell_field, has_cell_field2
    character(len=256) :: field_label, field_label2

    nvt = size(points, 2)
    if (present(elem_indices)) then
      nel = size(elem_indices)
    else
      nel = size(connectivity, 2)
    end if

    allocate(conn_flat(8*nel))
    allocate(offsets(nel))
    allocate(cellTypes(nel))

    has_cell_field = .false.
    has_cell_field2 = .false.
    if (present(cell_field)) then
      if (size(cell_field) /= nel) then
        write(*,*) 'Cell field size mismatch when writing ', trim(filename)
        deallocate(conn_flat, offsets, cellTypes)
        return
      end if
      has_cell_field = .true.
      if (present(field_name)) then
        field_label = trim(field_name)
      else
        field_label = 'cell_field'
      end if
    end if
    if (present(cell_field2)) then
      if (size(cell_field2) /= nel) then
        write(*,*) 'Second cell field size mismatch when writing ', trim(filename)
        deallocate(conn_flat, offsets, cellTypes)
        return
      end if
      has_cell_field2 = .true.
      if (present(field_name2)) then
        field_label2 = trim(field_name2)
      else
        field_label2 = 'cell_field_2'
      end if
    end if

    if (present(elem_indices)) then
      do i = 1, nel
        idx = elem_indices(i)
        conn_flat(8*(i-1)+1:8*i) = connectivity(:, idx) - 1
      end do
    else
      conn_flat = reshape(connectivity, [8*nel]) - 1
    end if
    offsets = [(8*i, i=1,nel)]
    cellTypes = 12

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Failed to open file ', trim(filename)
      deallocate(conn_flat, offsets, cellTypes)
      return
    end if

    write(unit,'(A)') '<?xml version="1.0"?>'
    write(unit,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(unit,'(A)') '  <UnstructuredGrid>'
    write(unit,'(A," NumberOfPoints=""",I0,""" NumberOfCells=""",I0,""">")') '    <Piece', nvt, nel
    write(unit,'(A)') '      <PointData/>'
    if (has_cell_field .or. has_cell_field2) then
      write(unit,'(A)') '      <CellData>'
      if (has_cell_field) then
        write(unit,'(A," type=""Float64"" Name=""",A,""" format=""ascii"">")') '        <DataArray', trim(field_label)
        do i = 1, nel
          write(unit,'(A,ES23.15)') '          ', cell_field(i)
        end do
        write(unit,'(A)') '        </DataArray>'
      end if
      if (has_cell_field2) then
        write(unit,'(A," type=""Float64"" Name=""",A,""" format=""ascii"">")') '        <DataArray', trim(field_label2)
        do i = 1, nel
          write(unit,'(A,ES23.15)') '          ', cell_field2(i)
        end do
        write(unit,'(A)') '        </DataArray>'
      end if
      write(unit,'(A)') '      </CellData>'
    else
      write(unit,'(A)') '      <CellData/>'
    end if
    write(unit,'(A)') '      <Points>'
    write(unit,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    do i = 1, nvt
      write(unit,'(A,3(1X,ES23.15))') '          ', points(1,i), points(2,i), points(3,i)
    end do
    write(unit,'(A)') '        </DataArray>'
    write(unit,'(A)') '      </Points>'
    write(unit,'(A)') '      <Cells>'
    write(unit,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
    write(unit,'(A,*(I0,1X))') '          ', conn_flat
    write(unit,'(A)') '        </DataArray>'
    write(unit,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
    write(unit,'(A,*(I0,1X))') '          ', offsets
    write(unit,'(A)') '        </DataArray>'
    write(unit,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
    write(unit,'(A,*(I0,1X))') '          ', cellTypes
    write(unit,'(A)') '        </DataArray>'
    write(unit,'(A)') '      </Cells>'
    write(unit,'(A)') '    </Piece>'
    write(unit,'(A)') '  </UnstructuredGrid>'
    write(unit,'(A)') '</VTKFile>'
    close(unit)

    deallocate(conn_flat, offsets, cellTypes)
  end subroutine write_hex_mesh_vtu

  subroutine write_pvtu_reference(path, nPieces, prefix, suffix, cell_field, cell_fields)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nPieces
    character(len=*), intent(in) :: prefix, suffix
    character(len=*), intent(in), optional :: cell_field
    character(len=*), intent(in), optional :: cell_fields(:)

    integer :: unit, ios, i, f, n_fields
    logical :: has_field, use_array_fields
    character(len=512) :: piece_file
    character(len=256), allocatable :: field_names(:)
    character(len=512) :: prefix_base
    integer :: lenf, slashpos, ic

    has_field = .false.
    use_array_fields = .false.
    n_fields = 0
    if (present(cell_fields)) then
      n_fields = size(cell_fields)
      if (n_fields > 0) then
        has_field = .true.
        use_array_fields = .true.
        allocate(field_names(n_fields))
        do i = 1, n_fields
          field_names(i) = trim(cell_fields(i))
        end do
      end if
    else if (present(cell_field)) then
      if (len_trim(cell_field) > 0) then
        has_field = .true.
        n_fields = 1
        allocate(field_names(1))
        field_names(1) = trim(cell_field)
      end if
    end if

    if (has_field) then
      if (n_fields == 0) has_field = .false.
    else
      n_fields = 0
    end if

    open(newunit=unit, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Failed to open PVTU file ', trim(path)
      return
    end if

    prefix_base = trim(prefix)
    lenf = len_trim(prefix_base)
    slashpos = 0
    do ic = lenf, 1, -1
      if (prefix_base(ic:ic) == '/' .or. prefix_base(ic:ic) == '\') then
        slashpos = ic
        exit
      end if
    end do
    if (slashpos > 0) prefix_base = prefix_base(slashpos+1:)

    write(unit,'(A)') '<?xml version="1.0"?>'
    write(unit,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(unit,'(A)') '  <PUnstructuredGrid GhostLevel="0">'
    if (has_field) then
      write(unit,'(A)') '    <PCellData>'
      do f = 1, n_fields
        if (use_array_fields) then
          write(unit,'(A," type=""Float64"" Name=""",A,""" NumberOfComponents=""1""/>")') '      <PDataArray', &
               trim(field_names(f))
        else
          write(unit,'(A," type=""Int32"" Name=""",A,"""/>")') '      <PDataArray', trim(field_names(f))
        end if
      end do
      write(unit,'(A)') '    </PCellData>'
    else
      write(unit,'(A)') '    <PCellData/>'
    end if
    write(unit,'(A)') '    <PPoints>'
    write(unit,'(A)') '      <PDataArray type="Float64" NumberOfComponents="3"/>'
    write(unit,'(A)') '    </PPoints>'
    write(unit,'(A)') '    <PCells>'
    write(unit,'(A)') '      <PDataArray type="Int32" Name="connectivity"/>'
    write(unit,'(A)') '      <PDataArray type="Int32" Name="offsets"/>'
    write(unit,'(A)') '      <PDataArray type="UInt8" Name="types"/>'
    write(unit,'(A)') '    </PCells>'
    do i = 0, nPieces-1
      write(piece_file,'(A,I4.4,A)') trim(prefix_base), i, trim(suffix)
      write(unit,'(A," Source=""",A,"""/>")') '    <Piece', trim(piece_file)
    end do
    write(unit,'(A)') '  </PUnstructuredGrid>'
    write(unit,'(A)') '</VTKFile>'
    close(unit)
    if (allocated(field_names)) deallocate(field_names)
  end subroutine write_pvtu_reference

end module hex_io
