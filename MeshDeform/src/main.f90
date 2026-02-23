program meshdeform_main
  use def_mod, only: hex_mesh_type, deformation_context, initialize_deformation_context, &
                     finalize_deformation_context, summarize_mesh, rk, assign_vertex_constraints, &
                     report_constraint_summary, vertex_constraint_lists, build_constraint_lists, &
                     release_constraint_lists, surface_mesh, load_surface_mesh, report_surface_summary, &
                     identify_boundary_nodes, build_edge_to_element, report_boundary_assignment_status, &
                     report_distance_summary, apply_edge_deformation, write_deformed_vtu
  use inout_mod, only: load_tri_mesh
  use preprocessor_config_mod, only: get_monitor_threshold_default
  use setupe3dfile_reader, only: MeshConfig, ProcessParameters, initialize_mesh_config, load_mesh_config, &
                                 log_mesh_config, initialize_process_parameters, load_process_parameters, &
                                 log_process_inflows
  implicit none

  character(len=1024) :: working_folder
  character(len=1024) :: input_mesh_path
  character(len=1024) :: setup_file_path
  character(len=1024) :: surface_file_path
  logical :: verbose
  type(hex_mesh_type) :: mesh
  type(deformation_context) :: context
  type(surface_mesh) :: surface
  type(MeshConfig) :: mesh_config
  type(ProcessParameters) :: process_params
  integer, allocatable :: vertex_constraints(:)
  type(vertex_constraint_lists) :: constraint_lists
  logical, allocatable :: boundary_mask(:)
  real(rk) :: monitor_threshold
  logical :: has_setup

  call print_banner('MESHDEFORM START')
  call parse_command_line(working_folder, input_mesh_path, verbose)
  monitor_threshold = real(get_monitor_threshold_default(), rk)

  if (len_trim(working_folder) == 0) then
     call print_usage('Working folder (-f) is required.', .true.)
  end if

  if (len_trim(input_mesh_path) == 0) then
     if (len_trim(working_folder) == 0) then
        call print_usage('Either an input mesh or folder must be provided.', .true.)
     end if
     input_mesh_path = trim(working_folder) // '/Filtered.tri'
  end if

  setup_file_path = trim(working_folder) // '/setup.e3d'
  inquire(file=trim(setup_file_path), exist=has_setup)
  if (.not.has_setup) then
     write(*, '(A)') 'Setup file not found: ' // trim(setup_file_path)
     stop 1
  end if

  call initialize_mesh_config(mesh_config)
  call initialize_process_parameters(process_params)
  call load_mesh_config(trim(setup_file_path), mesh_config)
  if (.not.mesh_config%loaded) then
     write(*, '(A)') 'Failed to load mesh configuration from: ' // trim(setup_file_path)
     stop 1
  end if
  call load_process_parameters(trim(setup_file_path), process_params)
  if (.not.process_params%loaded) then
     write(*, '(A)') 'Failed to load process parameters from: ' // trim(setup_file_path)
     stop 1
  end if

  surface_file_path = trim(working_folder) // '/surface.off'
  call load_surface_mesh(surface, trim(surface_file_path))
  call report_surface_summary(surface)

  call initialize_deformation_context(context)
  call load_tri_mesh(trim(input_mesh_path), mesh)
  if (mesh%nel <= 0 .or. mesh%nvt <= 0) then
     write(*, '(A)') 'Mesh file did not contain any elements; aborting.'
     call finalize_deformation_context(context, mesh)
     stop 1
  end if
  call identify_boundary_nodes(mesh, boundary_mask)

  if (verbose) then
     write(*, '(A,F8.3)') 'Monitor threshold default: ', monitor_threshold
     call log_mesh_config(mesh_config)
     call log_process_inflows(process_params)
  end if

  if (.not.verbose) then
     write(*, '(A)') 'Loaded mesh metadata from setup.e3d.'
  end if

  call assign_vertex_constraints(mesh, mesh_config, vertex_constraints)
  call build_constraint_lists(vertex_constraints, constraint_lists, boundary_mask)
  call report_constraint_summary(vertex_constraints, boundary_mask)
  call report_boundary_assignment_status(vertex_constraints, boundary_mask)

  call build_edge_to_element(mesh)
  call report_distance_summary(surface, mesh)
  call apply_edge_deformation(mesh, mesh_config, surface, 1)
  call write_deformed_vtu(mesh, surface, trim(working_folder) // '/MeshDeformResult.vtu')

  call summarize_mesh(mesh, trim(input_mesh_path))
  call finalize_deformation_context(context, mesh)
  if (allocated(vertex_constraints)) deallocate(vertex_constraints)
  if (allocated(boundary_mask)) deallocate(boundary_mask)
  call release_constraint_lists(constraint_lists)
  call surface%clear()
  call print_banner('MESHDEFORM END')

contains

  subroutine parse_command_line(working_folder, input_mesh, verbose)
    character(len=*), intent(out) :: working_folder
    character(len=*), intent(out) :: input_mesh
    logical, intent(out) :: verbose
    integer :: argc, idx
    character(len=1024) :: arg

    working_folder = ''
    input_mesh = ''
    verbose = .false.
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
       case ('-i', '--input')
          if (idx == argc) then
             call print_usage('Missing filename after input flag.', .true.)
          end if
          idx = idx + 1
          call get_command_argument(idx, arg)
          input_mesh = trim(arg)
       case ('-v', '--verbose')
          verbose = .true.
       case ('-h', '--help')
          call print_usage('', .false.)
       case default
          call print_usage('Unknown argument: ' // trim(arg), .true.)
       end select
       idx = idx + 1
    end do
  end subroutine parse_command_line

  subroutine print_usage(message, is_error)
    character(len=*), intent(in) :: message
    logical, intent(in) :: is_error

    if (len_trim(message) > 0) then
       write(*, '(A)') trim(message)
    end if
    write(*, '(A)') 'Usage: meshdeform -f <folder> [-i <mesh>] [-v]'
    write(*, '(A)') '       meshdeform --folder <folder> [--input <mesh>] [--verbose]'
    write(*, '(A)') 'The working folder must contain setup.e3d (mesh metadata) and typically Filtered.tri.'
    if (is_error) then
       stop 1
    else
       stop 0
    end if
  end subroutine print_usage

  subroutine print_banner(label)
    character(len=*), intent(in) :: label
    integer, parameter :: total_width = 96
    integer :: padding, left_pad, right_pad
    character(len=total_width) :: line

    padding = total_width - len_trim(label) - 2
    if (padding < 0) padding = 0
    left_pad = padding / 2
    right_pad = padding - left_pad
    line = repeat('=', left_pad) // ' ' // trim(label) // ' ' // repeat('=', right_pad)
    write(*, '(A)') line
  end subroutine print_banner

end program meshdeform_main
