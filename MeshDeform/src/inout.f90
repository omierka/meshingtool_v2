module inout_mod
  use iso_fortran_env, only: error_unit
  use def_mod, only: hex_mesh_type, rk
  implicit none
  private

  public :: load_tri_mesh

contains

  subroutine load_tri_mesh(filename, mesh)
    character(len=*), intent(in) :: filename
    type(hex_mesh_type), intent(inout) :: mesh

    integer :: unit, ios, iel, ivt
    character(len=256) :: line
    logical :: exists

    call mesh%clear()
    inquire(file=trim(filename), exist=exists)
    if (.not.exists) then
       write(error_unit, '(A)') 'Mesh file not found: ' // trim(filename)
       return
    end if

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(error_unit, '(A,I0)') 'Unable to open mesh file: ' // trim(filename) // ' IOSTAT=', ios
       return
    end if

    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) then
       call abort_read('Failed to read first header line from ', filename, unit, mesh)
       return
    end if
    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) then
       call abort_read('Failed to read second header line from ', filename, unit, mesh)
       return
    end if

    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) then
          call abort_read('Failed to read mesh counts from ', filename, unit, mesh)
          return
       end if
       if (len_trim(line) == 0) cycle
       if (index(adjustl(line), 'NEL') == 1) cycle
       exit
    end do

    read(line, *, iostat=ios) mesh%nel, mesh%nvt, mesh%nbct, mesh%nve, mesh%nee, mesh%nae
    if (ios /= 0) then
       call abort_read('Failed to read mesh counts from ', filename, unit, mesh)
       return
    end if

    call expect_section(unit, 'DCORVG', filename, mesh)
    allocate(mesh%dcorvg(3, mesh%nvt))
    do ivt = 1, mesh%nvt
       read(unit, *, iostat=ios) mesh%dcorvg(1, ivt), mesh%dcorvg(2, ivt), mesh%dcorvg(3, ivt)
       if (ios /= 0) then
          call abort_read('Failed to read DCORVG data from ', filename, unit, mesh)
          return
       end if
    end do

    call expect_section(unit, 'KVERT', filename, mesh)
    allocate(mesh%kvert(mesh%nve, mesh%nel))
    do iel = 1, mesh%nel
       read(unit, *, iostat=ios) mesh%kvert(:, iel)
       if (ios /= 0) then
          call abort_read('Failed to read KVERT connectivity from ', filename, unit, mesh)
          return
       end if
    end do

    call expect_section(unit, 'KNPR', filename, mesh)
    allocate(mesh%knpr(mesh%nvt))
    do ivt = 1, mesh%nvt
       read(unit, *, iostat=ios) mesh%knpr(ivt)
       if (ios /= 0) then
          call abort_read('Failed to read KNPR data from ', filename, unit, mesh)
          return
       end if
    end do

    close(unit)
  end subroutine load_tri_mesh

  subroutine expect_section(unit, expected, filename, mesh)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: expected
    character(len=*), intent(in) :: filename
    type(hex_mesh_type), intent(inout) :: mesh
    integer :: ios
    character(len=256) :: tag

    read(unit, '(A)', iostat=ios) tag
    if (ios /= 0) then
       call abort_read('Failed to read section marker from ', filename, unit, mesh)
       return
    end if
    if (.not.section_matches(tag, expected)) then
       write(error_unit, '(A)') 'Warning: Expected section [' // trim(expected) // '] but found [' // trim(tag) // ']'
    end if
  end subroutine expect_section

  logical function section_matches(tag, expected)
    character(len=*), intent(in) :: tag
    character(len=*), intent(in) :: expected
    character(len=:), allocatable :: normalized_tag, normalized_expected

    normalized_tag = adjustl(tag)
    normalized_expected = adjustl(expected)
    section_matches = (len_trim(normalized_tag) >= len_trim(normalized_expected)) .and. &
         (normalized_tag(1:len_trim(normalized_expected)) == normalized_expected(1:len_trim(normalized_expected)))
  end function section_matches

  subroutine abort_read(message, filename, unit, mesh)
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: filename
    integer, intent(in) :: unit
    type(hex_mesh_type), intent(inout) :: mesh

    write(error_unit, '(A)') trim(message) // trim(filename)
    close(unit)
    call mesh%clear()
  end subroutine abort_read

end module inout_mod
