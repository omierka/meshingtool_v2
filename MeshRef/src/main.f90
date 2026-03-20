program meshref_main
  use def_mod, only: initialize_mesh_database, finalize_mesh_database, &
                     load_target_mesh_file, build_target_element_span, build_face_adjacency, mesh_refinement, &
                     random_element_marking, vertice_marking, enforce_refinement_levels, load_monitor_file, &
                     report_refinement_distribution, mark_elements_by_threshold, apply_inflow_refinement_boost
  use inout_mod, only: write_patch_group_vtu, build_clean_output_path, &
                       build_refined_clean_output_path, build_level_output_path, write_refined_clean_tri
  use cleanup_mod, only: clear_intra_patches, clear_inter_patches, write_refined_clean_vtu, &
                         update_coord_tolerance_from_patches
  use var_mod, only: element_patches, clean_element_patches, hex_mesh, initialize_mesh_levels, &
                     bind_refined_mesh, default_refinement_depth, rk, Monitor_threshold
  use preprocessor_config_mod, only: get_monitor_threshold_default
  implicit none

  character(len=1024) :: working_folder
  character(len=1024) :: input_mesh_path
  character(len=1024) :: output_vtu_path
  character(len=1024) :: level_output_vtu_path
  character(len=1024) :: clean_output_vtu_path
  character(len=1024) :: refined_clean_output_vtu_path
  character(len=1024) :: monitor_file_path
  character(len=1024) :: setup_file_path
  character(len=1024) :: label
  integer :: random_percentage
  integer :: recursion_depth
  integer :: level
  real(rk) :: level_threshold
  logical :: inflow_levels_changed

  call print_banner('MESHREF START')
  recursion_depth = default_refinement_depth
  Monitor_threshold = real(get_monitor_threshold_default(), rk)
  call parse_command_line(working_folder, random_percentage, recursion_depth)
  if (random_percentage >= 0) then
     write(*, '(A)') 'Random refinement (-r/--random-refinement) is not supported when inflow boosting is enabled.'
     stop 1
  end if
  input_mesh_path = trim(working_folder) // '/Coarse_meshDir/Mesh.tri'
  monitor_file_path = trim(working_folder) // '/area.txt'
  setup_file_path = trim(working_folder) // '/setup.e3d'
  output_vtu_path = trim(working_folder) // '/RefinedCleanMesh.vtu'

  call initialize_mesh_database()
  call initialize_mesh_levels(recursion_depth)
  call load_target_mesh_file(trim(input_mesh_path))
  call build_target_element_span(hex_mesh(0))
  call build_face_adjacency(hex_mesh(0))
  call load_monitor_file(hex_mesh(0), trim(monitor_file_path))
  if (random_percentage >= 0) call random_element_marking(random_percentage)
  level_threshold = Monitor_threshold
  do level = recursion_depth, 1, -1
     if (random_percentage < 0) then
        call mark_elements_by_threshold(hex_mesh(0), level, level_threshold)
     end if
     call enforce_refinement_levels(hex_mesh(0), level)
     if (random_percentage < 0 .and. level > 1) level_threshold = level_threshold / 3.0_rk
  end do
  call apply_inflow_refinement_boost(hex_mesh(0), trim(setup_file_path), inflow_levels_changed)
  if (inflow_levels_changed) then
     do level = recursion_depth, 1, -1
        call enforce_refinement_levels(hex_mesh(0), level)
     end do
  end if
  write(label, '(A,I0,A)') 'No Of Elements with Refinement depth [0..', recursion_depth, '] (smoothed)'
  call report_refinement_distribution(hex_mesh(0), trim(label))

  do level = 0, recursion_depth - 1

     if (.not.allocated(hex_mesh(level)%kelementspan)) call build_target_element_span(hex_mesh(level))
     call vertice_marking(hex_mesh(level), level+1)

     if (level == 0) then
        call write_refined_clean_vtu(hex_mesh(level), trim('target.vtu'))
     end if

     call build_level_output_path(trim(output_vtu_path), level, level_output_vtu_path)
     call mesh_refinement(hex_mesh(level))
     if (.not.allocated(element_patches)) then
        write(*, '(A,I0)') 'No element patches present at refinement level ', level
        cycle
     end if
     call write_patch_group_vtu(element_patches, trim(level_output_vtu_path))

     call update_coord_tolerance_from_patches(element_patches)
     call clear_intra_patches(element_patches)
     if (.not.allocated(clean_element_patches)) cycle

     call build_clean_output_path(trim(level_output_vtu_path), clean_output_vtu_path)
     call write_patch_group_vtu(clean_element_patches, trim(clean_output_vtu_path))

     call bind_refined_mesh(level + 1)
     call clear_inter_patches(hex_mesh(level))
     call build_target_element_span(hex_mesh(level + 1))
     call build_refined_clean_output_path(trim(level_output_vtu_path), refined_clean_output_vtu_path)
     call write_refined_clean_vtu(hex_mesh(level + 1), trim(refined_clean_output_vtu_path))

  end do

  write(label, '(A,I0,A)') 'Final No Of Elements with Refinement depth [0..', recursion_depth, ']'
  call report_refinement_distribution(hex_mesh(recursion_depth), trim(label))

  call write_refined_clean_tri(hex_mesh(recursion_depth), hex_mesh(0), trim(working_folder))

  call finalize_mesh_database()
  call print_banner('MESHREF END')

contains

  subroutine print_banner(label)
    character(len=*), intent(in) :: label
    integer, parameter :: total_width = 110
    integer :: padding, left_width, right_width
    character(len=total_width) :: line

    padding = total_width - len_trim(label) - 2
    if (padding < 0) padding = 0
    left_width = padding / 2
    right_width = padding - left_width
    line = repeat('=', left_width) // ' ' // trim(label) // ' ' // repeat('=', right_width)
    write(*, '(A)') line
  end subroutine print_banner

  subroutine parse_command_line(working_folder, random_percentage, recursion_depth)
    use var_mod, only: rk, Monitor_threshold, default_refinement_depth
    character(len=*), intent(out) :: working_folder
    integer, intent(out) :: random_percentage
    integer, intent(out) :: recursion_depth
    integer :: argc, idx, ios, value
    character(len=1024) :: arg

    working_folder = ''
    random_percentage = -1
    recursion_depth = default_refinement_depth
    argc = command_argument_count()
    idx = 1

    do while (idx <= argc)
       call get_command_argument(idx, arg)
      select case (trim(arg))
       case ('-f', '--folder')
         if (idx == argc) then
             call print_usage('Missing folder path after folder flag.', .true.)
          end if
          idx = idx + 1
          call get_command_argument(idx, arg)
          working_folder = trim(arg)
       case ('-r', '--random-refinement')
          if (idx == argc) then
             call print_usage('Missing percentage after random refinement flag.', .true.)
          end if
          idx = idx + 1
          call get_command_argument(idx, arg)
          read(arg, *, iostat=ios) value
          if (ios /= 0 .or. value < 0 .or. value > 100) then
             call print_usage('Random refinement value must be between 0 and 100.', .true.)
          end if
          random_percentage = value
       case ('-d', '--depth', '--recursion-depth')
          if (idx == argc) then
             call print_usage('Missing value after depth flag.', .true.)
          end if
          idx = idx + 1
          call get_command_argument(idx, arg)
          read(arg, *, iostat=ios) value
          if (ios /= 0 .or. value < 1 .or. value > 3) then
             call print_usage('Depth must be an integer between 1 and 3.', .true.)
          end if
          recursion_depth = value
       case ('-h', '--help')
          call print_usage('', .false.)
       case default
          call print_usage('Unknown argument: ' // trim(arg), .true.)
       end select
       idx = idx + 1
    end do

    if (len_trim(working_folder) == 0) then
       call print_usage('Working folder (-f) is required.', .true.)
    end if
  end subroutine parse_command_line

  subroutine print_usage(message, is_error)
    character(len=*), intent(in) :: message
    logical, intent(in) :: is_error

    if (len_trim(message) > 0) then
       write(*, '(A)') trim(message)
    end if
    write(*, '(A)') 'Usage: meshref -f <folder> [options]'
    write(*, '(A)') '       meshref --folder <folder> [options]'
    write(*, '(A)') 'Options:'
    write(*, '(A)') '  -r, --random-refinement <0-100>  (disabled; area.txt + inflow boost required).'
    write(*, '(A)') '  -d, --depth <1-3>                Recursion depth (default 2).'
    write(*, '(A)') '  -h, --help                       Show this help text.'
    if (is_error) then
       stop 1
    else
       stop 0
    end if
  end subroutine print_usage

end program meshref_main
