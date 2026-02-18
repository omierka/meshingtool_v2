module preprocessor_config_mod
  use iso_fortran_env, only: real64, error_unit
  implicit none
  private

  character(len=*), parameter :: default_config_filename = "preprocessor.cfg"

  logical :: config_loaded = .false.

  real(real64) :: cfg_cylinder_boundary_tolerance = 0.90_real64
  real(real64) :: cfg_box_boundary_tolerance = 0.50_real64
  real(real64) :: cfg_tolerance_floor = 1.0e-12_real64
  real(real64) :: cfg_inflow_cos_threshold = 0.99_real64
  real(real64) :: cfg_monitor_threshold_default = 1.5_real64
  real(real64) :: cfg_triangle_tree_tol = 1.0e-12_real64
  real(real64) :: cfg_clip_eps = 1.0e-12_real64

  public :: get_cylinder_boundary_tolerance_factor
  public :: get_box_boundary_tolerance_factor
  public :: get_tolerance_floor
  public :: get_inflow_orientation_threshold
  public :: get_monitor_threshold_default
  public :: get_triangle_tree_tolerance
  public :: get_monitor_clip_epsilon

contains

  subroutine ensure_config_loaded()
    if (config_loaded) return
    call load_config()
    config_loaded = .true.
  end subroutine ensure_config_loaded

  subroutine load_config()
    character(len=:), allocatable :: path
    character(len=1024) :: raw_line
    character(len=256) :: section
    character(len=256) :: key
    character(len=1024) :: value
    integer :: unit, ios, eq_pos
    character(len=:), allocatable :: trimmed

    call resolve_config_path(path)
    open(newunit=unit, file=trim(path), status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write(error_unit, '(A,1X,A)') "Warning: unable to open preprocessor config:", trim(path)
      return
    end if

    section = ""
    do
      read(unit, '(A)', iostat=ios) raw_line
      if (ios /= 0) exit
      trimmed = adjustl(raw_line)
      if (len_trim(trimmed) == 0) cycle
      if (trimmed(1:1) == "#" .or. trimmed(1:1) == ";") cycle
      if (trimmed(1:1) == "[") then
        call parse_section_name(trimmed, section)
        cycle
      end if
      eq_pos = index(trimmed, "=")
      if (eq_pos <= 1) cycle
      key = lowercase(trim(adjustl(trimmed(:eq_pos - 1))))
      value = trim(adjustl(trimmed(eq_pos + 1:)))
      call assign_value(trim(lowercase(section)), trim(key), value)
    end do
    close(unit)
  end subroutine load_config

  subroutine resolve_config_path(path)
    character(len=:), allocatable, intent(out) :: path
    character(len=1024) :: buffer
    integer :: length, status

    call get_environment_variable("PREPROCESSOR_CONFIG", buffer, length=length, status=status)
    if (status == 0 .and. length > 0) then
      path = buffer(:length)
    else
      path = default_config_filename
    end if
  end subroutine resolve_config_path

  subroutine assign_value(section, key, raw_value)
    character(len=*), intent(in) :: section
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: raw_value

    select case (trim(section))
    case ("3dmeshcleaner")
      select case (trim(key))
      case ("cylinder_boundary_tolerance_factor")
        call assign_real(raw_value, cfg_cylinder_boundary_tolerance)
      case ("box_boundary_tolerance_factor")
        call assign_real(raw_value, cfg_box_boundary_tolerance)
      case ("tolerance_floor")
        call assign_real(raw_value, cfg_tolerance_floor)
      case ("inflowcos_threshold_orientation_")
        call assign_real(raw_value, cfg_inflow_cos_threshold)
      end select
    case ("3dmeshref")
      if (trim(key) == "monitor_threshold_default") then
        call assign_real(raw_value, cfg_monitor_threshold_default)
      end if
    case ("3dmonitorgen")
      select case (trim(key))
      case ("triangle_tree_tol")
        call assign_real(raw_value, cfg_triangle_tree_tol)
      case ("clip_eps")
        call assign_real(raw_value, cfg_clip_eps)
      end select
    end select
  end subroutine assign_value

  subroutine assign_real(raw_value, target)
    character(len=*), intent(in) :: raw_value
    real(real64), intent(inout) :: target
    real(real64) :: parsed
    integer :: ios

    read(raw_value, *, iostat=ios) parsed
    if (ios == 0) then
      target = parsed
    end if
  end subroutine assign_real

  subroutine parse_section_name(line, section)
    character(len=*), intent(in) :: line
    character(len=*), intent(inout) :: section
    integer :: close_idx

    close_idx = index(line, "]")
    if (close_idx > 2) then
      section = trim(line(2:close_idx - 1))
    else
      section = ""
    end if
  end subroutine parse_section_name

  pure function lowercase(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: out
    integer :: idx, code

    do idx = 1, len(str)
      code = iachar(str(idx:idx))
      if (code >= iachar("A") .and. code <= iachar("Z")) then
        out(idx:idx) = achar(code + 32)
      else
        out(idx:idx) = str(idx:idx)
      end if
    end do
  end function lowercase

  real(real64) function get_cylinder_boundary_tolerance_factor()
    call ensure_config_loaded()
    get_cylinder_boundary_tolerance_factor = cfg_cylinder_boundary_tolerance
  end function get_cylinder_boundary_tolerance_factor

  real(real64) function get_box_boundary_tolerance_factor()
    call ensure_config_loaded()
    get_box_boundary_tolerance_factor = cfg_box_boundary_tolerance
  end function get_box_boundary_tolerance_factor

  real(real64) function get_tolerance_floor()
    call ensure_config_loaded()
    get_tolerance_floor = cfg_tolerance_floor
  end function get_tolerance_floor

  real(real64) function get_inflow_orientation_threshold()
    call ensure_config_loaded()
    get_inflow_orientation_threshold = cfg_inflow_cos_threshold
  end function get_inflow_orientation_threshold

  real(real64) function get_monitor_threshold_default()
    call ensure_config_loaded()
    get_monitor_threshold_default = cfg_monitor_threshold_default
  end function get_monitor_threshold_default

  real(real64) function get_triangle_tree_tolerance()
    call ensure_config_loaded()
    get_triangle_tree_tolerance = cfg_triangle_tree_tol
  end function get_triangle_tree_tolerance

  real(real64) function get_monitor_clip_epsilon()
    call ensure_config_loaded()
    get_monitor_clip_epsilon = cfg_clip_eps
  end function get_monitor_clip_epsilon

end module preprocessor_config_mod
