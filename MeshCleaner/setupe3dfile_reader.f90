module setupe3dfile_reader
    use iso_c_binding, only: c_double
    implicit none
    private

    type, public :: BoxMeshConfig
        real(c_double) :: geometry_start(3) = 0.0_c_double
        real(c_double) :: geometry_length(3) = 0.0_c_double
        real(c_double) :: spacing(3) = 0.0_c_double
        logical :: has_geometry_start = .false.
        logical :: has_geometry_length = .false.
        logical :: has_spacing = .false.
    end type BoxMeshConfig

    type, public :: HollowCylinderMeshConfig
        real(c_double) :: barrel_diameter = 0.0_c_double
        real(c_double) :: inner_diameter = 0.0_c_double
        real(c_double) :: barrel_length = 0.0_c_double
        real(c_double) :: axial_start = 0.0_c_double
        real(c_double) :: spacing(3) = 0.0_c_double
        logical :: has_barrel_diameter = .false.
        logical :: has_inner_diameter = .false.
        logical :: has_barrel_length = .false.
        logical :: has_axial_start = .false.
        logical :: has_spacing = .false.
    end type HollowCylinderMeshConfig

    type, public :: MeshConfig
        character(len=64) :: mesh_type = "UNSPECIFIED"
        logical :: loaded = .false.
        character(len=512) :: source_file = ""
        type(BoxMeshConfig) :: box
        type(HollowCylinderMeshConfig) :: cylinder
    end type MeshConfig

    type, public :: ProcessInflow
        integer :: type_id = 0
        character(len=32) :: type_label = ""
        real(c_double) :: inner_radius = 0.0_c_double
        real(c_double) :: outer_radius = 0.0_c_double
        real(c_double) :: center(3) = 0.0_c_double
        real(c_double) :: normal(3) = 0.0_c_double
        real(c_double) :: midpoint_a(3) = 0.0_c_double
        real(c_double) :: midpoint_b(3) = 0.0_c_double
        logical :: has_type = .false.
        logical :: has_inner_radius = .false.
        logical :: has_outer_radius = .false.
        logical :: has_center = .false.
        logical :: has_normal = .false.
        logical :: has_midpoint_a = .false.
        logical :: has_midpoint_b = .false.
    end type ProcessInflow

    type, public :: ProcessParameters
        integer :: nOfInflows = 0
        character(len=512) :: source_file = ""
        logical :: loaded = .false.
        type(ProcessInflow), allocatable :: inflows(:)
    end type ProcessParameters

    public :: initialize_mesh_config, load_mesh_config, log_mesh_config
    public :: initialize_process_parameters, load_process_parameters, log_process_inflows

contains
    subroutine initialize_mesh_config(config)
        type(MeshConfig), intent(out) :: config

        config%mesh_type = "UNSPECIFIED"
        config%loaded = .false.
        config%source_file = ""
        config%box = BoxMeshConfig()
        config%cylinder = HollowCylinderMeshConfig()
    end subroutine initialize_mesh_config

    subroutine load_mesh_config(filename, config)
        character(len=*), intent(in) :: filename
        type(MeshConfig), intent(inout) :: config

        logical :: exists
        integer :: unit, ios, eq_pos, line_len
        character(len=512) :: raw_line, line, section, key, value

        inquire(file=filename, exist=exists)
        if (.not. exists) then
            write(*, '(A)') "Mesh configuration file not found: "//trim(filename)
            return
        end if

        call initialize_mesh_config(config)
        config%source_file = trim(filename)
        section = ""

        open(newunit=unit, file=filename, status="old", action="read")
        do
            read(unit, '(A)', iostat=ios) raw_line
            if (ios /= 0) exit
            line = adjustl(trim(raw_line))
            line_len = len_trim(line)
            if (line_len == 0) cycle
            if (line(1:1) == "!" .or. line(1:1) == "#") cycle
            if (line(1:1) == "[") then
                call parse_section_name(line, section)
                cycle
            end if
            eq_pos = index(line, "=")
            if (eq_pos <= 1 .or. eq_pos >= line_len) cycle
            key = trim(adjustl(line(:eq_pos - 1)))
            value = trim(adjustl(line(eq_pos + 1:)))
            call assign_config_value(section, key, value, config)
        end do
        close(unit)
        config%loaded = .true.
    end subroutine load_mesh_config

    subroutine log_mesh_config(config)
        type(MeshConfig), intent(in) :: config
        character(len=32) :: type_token

        if (.not. config%loaded) return

        write(*,'(A)') "Mesh configuration summary:"
        write(*,'(A,1X,A)') "  Source:", trim(config%source_file)
        write(*,'(A,1X,A)') "  HexMesher:", trim(config%mesh_type)
        type_token = lowercase(trim(config%mesh_type))
        select case (trim(type_token))
        case ("box")
            if (config%box%has_geometry_start) &
                write(*,'(A,3(1X,ES12.5))') "  geometryStart:", config%box%geometry_start
            if (config%box%has_geometry_length) &
                write(*,'(A,3(1X,ES12.5))') "  geometryLength:", config%box%geometry_length
            if (config%box%has_spacing) &
                write(*,'(A,3(1X,ES12.5))') "  sEl(x,y,z):", config%box%spacing
        case ("hollowcylinder")
            if (config%cylinder%has_barrel_diameter) &
                write(*,'(A,1X,ES12.5)') "  BarrelDiameter:", config%cylinder%barrel_diameter
            if (config%cylinder%has_inner_diameter) &
                write(*,'(A,1X,ES12.5)') "  InnerDiameter :", config%cylinder%inner_diameter
            if (config%cylinder%has_barrel_length) &
                write(*,'(A,1X,ES12.5)') "  BarrelLength  :", config%cylinder%barrel_length
            if (config%cylinder%has_axial_start) &
                write(*,'(A,1X,ES12.5)') "  AxialStartPos :", config%cylinder%axial_start
            if (config%cylinder%has_spacing) &
                write(*,'(A,3(1X,ES12.5))') "  sEl(tan,rad,ax):", config%cylinder%spacing
        end select
    end subroutine log_mesh_config

    subroutine assign_config_value(section, key, value, config)
        character(len=*), intent(in) :: section, key, value
        type(MeshConfig), intent(inout) :: config
        character(len=128) :: section_id, key_id
        logical :: handled

        section_id = lowercase(trim(section))
        key_id = lowercase(trim(key))
        handled = .false.

        select case (section_id)
        case ("e3dgeometrydata/preprocessing")
            handled = assign_simulation_key(key_id, value, config)
            if (.not. handled) handled = assign_geometry_key(key_id, value, config)
        case default
            handled = .false.
        end select
    end subroutine assign_config_value

    logical function assign_simulation_key(key_id, value, config)
        character(len=*), intent(in) :: key_id, value
        type(MeshConfig), intent(inout) :: config
        logical :: parsed

        assign_simulation_key = .false.
        select case (key_id)
        case ("hexmesher")
            call set_mesh_type(config, value)
            assign_simulation_key = .true.
        case ("sel_x")
            parsed = try_parse_real(value, config%box%spacing(1))
            if (parsed) then
                config%box%has_spacing = .true.
                assign_simulation_key = .true.
            end if
        case ("sel_y")
            parsed = try_parse_real(value, config%box%spacing(2))
            if (parsed) then
                config%box%has_spacing = .true.
                assign_simulation_key = .true.
            end if
        case ("sel_z")
            parsed = try_parse_real(value, config%box%spacing(3))
            if (parsed) then
                config%box%has_spacing = .true.
                assign_simulation_key = .true.
            end if
        case ("sel_tangential")
            parsed = try_parse_real(value, config%cylinder%spacing(1))
            if (parsed) then
                config%cylinder%has_spacing = .true.
                assign_simulation_key = .true.
            end if
        case ("sel_radial")
            parsed = try_parse_real(value, config%cylinder%spacing(2))
            if (parsed) then
                config%cylinder%has_spacing = .true.
                assign_simulation_key = .true.
            end if
        case ("sel_axial")
            parsed = try_parse_real(value, config%cylinder%spacing(3))
            if (parsed) then
                config%cylinder%has_spacing = .true.
                assign_simulation_key = .true.
            end if
        case default
            assign_simulation_key = .false.
        end select
    end function assign_simulation_key

    logical function assign_geometry_key(key_id, value, config)
        character(len=*), intent(in) :: key_id, value
        type(MeshConfig), intent(inout) :: config
        real(c_double) :: temp_vec(3)
        logical :: parsed

        assign_geometry_key = .false.
        select case (key_id)
        case ("geometrystart")
            parsed = try_parse_real_vector(value, temp_vec)
            if (parsed) then
                config%box%geometry_start = temp_vec
                config%box%has_geometry_start = .true.
                assign_geometry_key = .true.
            end if
        case ("geometrylength")
            parsed = try_parse_real_vector(value, temp_vec)
            if (parsed) then
                config%box%geometry_length = temp_vec
                config%box%has_geometry_length = .true.
                assign_geometry_key = .true.
            end if
        case ("barreldiameter")
            parsed = try_parse_real(value, config%cylinder%barrel_diameter)
            if (parsed) then
                config%cylinder%has_barrel_diameter = .true.
                assign_geometry_key = .true.
            end if
        case ("innerdiameter")
            parsed = try_parse_real(value, config%cylinder%inner_diameter)
            if (parsed) then
                config%cylinder%has_inner_diameter = .true.
                assign_geometry_key = .true.
            end if
        case ("barrellength")
            parsed = try_parse_real(value, config%cylinder%barrel_length)
            if (parsed) then
                config%cylinder%has_barrel_length = .true.
                assign_geometry_key = .true.
            end if
        case ("axialstartposition")
            parsed = try_parse_real(value, config%cylinder%axial_start)
            if (parsed) then
                config%cylinder%has_axial_start = .true.
                assign_geometry_key = .true.
            end if
        case default
            assign_geometry_key = .false.
        end select
    end function assign_geometry_key

    subroutine set_mesh_type(config, raw_value)
        type(MeshConfig), intent(inout) :: config
        character(len=*), intent(in) :: raw_value
        character(len=64) :: token

        token = lowercase(trim(raw_value))
        select case (trim(token))
        case ("box")
            config%mesh_type = "Box"
        case ("hollowcylinder")
            config%mesh_type = "HollowCylinder"
        case default
            config%mesh_type = trim(raw_value)
        end select
    end subroutine set_mesh_type

    subroutine initialize_process_parameters(process)
        type(ProcessParameters), intent(out) :: process

        process%nOfInflows = 0
        process%source_file = ""
        process%loaded = .false.
        if (allocated(process%inflows)) deallocate(process%inflows)
    end subroutine initialize_process_parameters

    subroutine load_process_parameters(filename, process)
        character(len=*), intent(in) :: filename
        type(ProcessParameters), intent(inout) :: process

        logical :: exists
        integer :: unit, ios, eq_pos, line_len
        character(len=512) :: raw_line, line, section, key, value

        inquire(file=filename, exist=exists)
        if (.not. exists) then
            write(*, '(A)') "Process parameter file not found: "//trim(filename)
            return
        end if

        call initialize_process_parameters(process)
        process%source_file = trim(filename)
        section = ""

        open(newunit=unit, file=filename, status="old", action="read")
        do
            read(unit, '(A)', iostat=ios) raw_line
            if (ios /= 0) exit
            line = adjustl(trim(raw_line))
            line_len = len_trim(line)
            if (line_len == 0) cycle
            if (line(1:1) == "!" .or. line(1:1) == "#") cycle
            if (line(1:1) == "[") then
                call parse_section_name(line, section)
                cycle
            end if
            eq_pos = index(line, "=")
            if (eq_pos <= 1 .or. eq_pos >= line_len) cycle
            key = trim(adjustl(line(:eq_pos - 1)))
            value = trim(adjustl(line(eq_pos + 1:)))
            call assign_process_value(section, key, value, process)
        end do
        close(unit)
        process%loaded = .true.
    end subroutine load_process_parameters

    subroutine log_process_inflows(process)
        type(ProcessParameters), intent(in) :: process
        integer :: idx

        if (.not. process%loaded) return

        write(*,'(A)') "Process inflow summary:"
        write(*,'(A,1X,A)') "  Source:", trim(process%source_file)
        write(*,'(A,I0)') "  nOfInflows:", process%nOfInflows
        if (.not. allocated(process%inflows)) return
        do idx = 1, min(process%nOfInflows, size(process%inflows))
            call log_single_inflow(idx, process%inflows(idx))
        end do
    end subroutine log_process_inflows

    subroutine log_single_inflow(idx, inflow)
        integer, intent(in) :: idx
        type(ProcessInflow), intent(in) :: inflow

        write(*,'(A,I0)') "    Inflow index:", idx
        if (inflow%has_type) then
            if (len_trim(inflow%type_label) > 0) then
                write(*,'(A,I0,2A)') "      type:", inflow%type_id, " ("//trim(inflow%type_label)//")"
            else
                write(*,'(A,I0)') "      type:", inflow%type_id
            end if
        else
            write(*,'(A)') "      type: (n/a)"
        end if
        if (inflow%has_inner_radius) then
            write(*,'(A,1X,ES12.5)') "      innerRadius:", inflow%inner_radius
        else
            write(*,'(A)') "      innerRadius: (n/a)"
        end if
        if (inflow%has_outer_radius) then
            write(*,'(A,1X,ES12.5)') "      outerRadius:", inflow%outer_radius
        else
            write(*,'(A)') "      outerRadius: (n/a)"
        end if
        if (inflow%has_center) then
            write(*,'(A,3(1X,ES12.5))') "      center:", inflow%center
        else
            write(*,'(A)') "      center: (n/a)"
        end if
        if (inflow%has_normal) then
            write(*,'(A,3(1X,ES12.5))') "      normal:", inflow%normal
        else
            write(*,'(A)') "      normal: (n/a)"
        end if
        if (inflow%has_midpoint_a) then
            write(*,'(A,3(1X,ES12.5))') "      midpointA:", inflow%midpoint_a
        else
            write(*,'(A)') "      midpointA: (n/a)"
        end if
        if (inflow%has_midpoint_b) then
            write(*,'(A,3(1X,ES12.5))') "      midpointB:", inflow%midpoint_b
        else
            write(*,'(A)') "      midpointB: (n/a)"
        end if
    end subroutine log_single_inflow

    subroutine assign_process_value(section, key, value, process)
        character(len=*), intent(in) :: section, key, value
        type(ProcessParameters), intent(inout) :: process
        character(len=128) :: section_id, key_id
        integer :: inflow_idx
        logical :: handled

        section_id = lowercase(trim(section))
        key_id = lowercase(trim(key))
        handled = .false.

        select case (section_id)
        case ("e3dprocessparameters")
            handled = assign_process_root_value(key_id, value, process)
        case default
            inflow_idx = parse_inflow_section_index(section_id)
            if (inflow_idx > 0) then
                handled = assign_inflow_key(inflow_idx, key_id, value, process)
            end if
        end select
    end subroutine assign_process_value

    logical function assign_process_root_value(key_id, value, process)
        character(len=*), intent(in) :: key_id, value
        type(ProcessParameters), intent(inout) :: process
        integer :: parsed_int
        logical :: parsed

        assign_process_root_value = .false.
        select case (key_id)
        case ("nofinflows")
            parsed = try_parse_integer(value, parsed_int)
            if (parsed .and. parsed_int >= 0) then
                call set_inflow_count(process, parsed_int)
                assign_process_root_value = .true.
            end if
        case default
            assign_process_root_value = .false.
        end select
    end function assign_process_root_value

    logical function assign_inflow_key(inflow_idx, key_id, value, process)
        integer, intent(in) :: inflow_idx
        character(len=*), intent(in) :: key_id, value
        type(ProcessParameters), intent(inout) :: process
        real(c_double) :: vec(3), scalar
        integer :: type_id
        character(len=32) :: label
        logical :: parsed

        assign_inflow_key = .false.
        if (.not. allocated(process%inflows)) return
        if (inflow_idx < 1 .or. inflow_idx > size(process%inflows)) return

        associate(inflow => process%inflows(inflow_idx))
            select case (key_id)
            case ("type")
                parsed = try_parse_inflow_type(value, type_id, label)
                if (parsed) then
                    inflow%type_id = type_id
                    inflow%type_label = label
                    inflow%has_type = .true.
                    assign_inflow_key = .true.
                end if
            case ("innerradius")
                parsed = try_parse_real(value, scalar)
                if (parsed) then
                    inflow%inner_radius = scalar
                    inflow%has_inner_radius = .true.
                    assign_inflow_key = .true.
                end if
            case ("outerradius")
                parsed = try_parse_real(value, scalar)
                if (parsed) then
                    inflow%outer_radius = scalar
                    inflow%has_outer_radius = .true.
                    assign_inflow_key = .true.
                end if
            case ("center")
                parsed = try_parse_real_vector(value, vec)
                if (parsed) then
                    inflow%center = vec
                    inflow%has_center = .true.
                    assign_inflow_key = .true.
                end if
            case ("normal")
                parsed = try_parse_real_vector(value, vec)
                if (parsed) then
                    inflow%normal = vec
                    inflow%has_normal = .true.
                    assign_inflow_key = .true.
                end if
            case ("midpointa")
                parsed = try_parse_real_vector(value, vec)
                if (parsed) then
                    inflow%midpoint_a = vec
                    inflow%has_midpoint_a = .true.
                    assign_inflow_key = .true.
                end if
            case ("midpointb")
                parsed = try_parse_real_vector(value, vec)
                if (parsed) then
                    inflow%midpoint_b = vec
                    inflow%has_midpoint_b = .true.
                    assign_inflow_key = .true.
                end if
            case default
                assign_inflow_key = .false.
            end select
        end associate
    end function assign_inflow_key

    subroutine set_inflow_count(process, count)
        type(ProcessParameters), intent(inout) :: process
        integer, intent(in) :: count

        if (allocated(process%inflows)) then
            deallocate(process%inflows)
        end if
        if (count > 0) then
            allocate(process%inflows(count))
        end if
        process%nOfInflows = count
    end subroutine set_inflow_count

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
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                out(idx:idx) = achar(code + 32)
            else
                out(idx:idx) = str(idx:idx)
            end if
        end do
    end function lowercase

    logical function try_parse_real(str, result_value)
        character(len=*), intent(in) :: str
        real(c_double), intent(out) :: result_value
        integer :: ios

        read(str, *, iostat=ios) result_value
        try_parse_real = (ios == 0)
    end function try_parse_real

    logical function try_parse_real_vector(str, values)
        character(len=*), intent(in) :: str
        real(c_double), intent(out) :: values(:)
        character(len=512) :: sanitized
        integer :: ios, idx

        sanitized = str
        do idx = 1, len(sanitized)
            if (sanitized(idx:idx) == ",") sanitized(idx:idx) = " "
        end do
        read(sanitized, *, iostat=ios) values
        try_parse_real_vector = (ios == 0)
    end function try_parse_real_vector

    logical function try_parse_integer(str, result_value)
        character(len=*), intent(in) :: str
        integer, intent(out) :: result_value
        integer :: ios

        read(str, *, iostat=ios) result_value
        try_parse_integer = (ios == 0)
    end function try_parse_integer

    logical function try_parse_inflow_type(str, type_id, label)
        character(len=*), intent(in) :: str
        integer, intent(out) :: type_id
        character(len=*), intent(out) :: label
        character(len=64) :: token
        logical :: parsed_int

        label = trim(str)
        token = lowercase(trim(str))
        select case (token)
        case ("rotatedparabola1")
            type_id = 1
            label = "ROTATEDPARABOLA1"
            try_parse_inflow_type = .true.
        case ("rotatedparabola2")
            type_id = 2
            label = "ROTATEDPARABOLA2"
            try_parse_inflow_type = .true.
        case ("flat")
            type_id = 3
            label = "FLAT"
            try_parse_inflow_type = .true.
        case ("curvedflat")
            type_id = 4
            label = "CURVEDFLAT"
            try_parse_inflow_type = .true.
        case ("rectangle")
            type_id = 5
            label = "RECTANGLE"
            try_parse_inflow_type = .true.
        case ("curvedrectangle")
            type_id = 6
            label = "CURVEDRECTANGLE"
            try_parse_inflow_type = .true.
        case default
            parsed_int = try_parse_integer(str, type_id)
            try_parse_inflow_type = parsed_int
        end select
        if (.not. try_parse_inflow_type) then
            type_id = 0
        end if
    end function try_parse_inflow_type

    integer function parse_inflow_section_index(section_id) result(idx)
        character(len=*), intent(in) :: section_id
        character(len=*), parameter :: prefix = "e3dprocessparameters/inflow_"
        integer :: ios

        idx = 0
        if (len_trim(section_id) < len(prefix)) return
        if (section_id(1:len(prefix)) /= prefix) return
        read(section_id(len(prefix)+1:), *, iostat=ios) idx
        if (ios /= 0) idx = 0
    end function parse_inflow_section_index

end module setupe3dfile_reader
