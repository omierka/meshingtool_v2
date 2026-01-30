module cgal_bindings
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_double, c_loc, c_associated
  implicit none
  private
  public :: cgal_initialize_triangle_tree
  public :: cgal_finalize_triangle_tree
  public :: cgal_have_triangle_tree
  public :: cgal_triangles_in_bbox

  logical, save :: tree_ready = .false.
  type(c_ptr), save :: tree_handle = c_null_ptr

  interface
    function cgal_create_triangle_tree(points, n_vertices, triangles, n_triangles) result(handle) &
         bind(C, name="cgal_create_triangle_tree")
      import :: c_ptr, c_int
      type(c_ptr) :: handle
      type(c_ptr), value :: points
      integer(c_int), value :: n_vertices
      type(c_ptr), value :: triangles
      integer(c_int), value :: n_triangles
    end function cgal_create_triangle_tree

    subroutine cgal_free_triangle_tree(handle) bind(C, name="cgal_free_triangle_tree")
      import :: c_ptr
      type(c_ptr), value :: handle
    end subroutine cgal_free_triangle_tree

    function cgal_count_triangles_in_bbox(handle, bbox_min, bbox_max) result(count) &
         bind(C, name="cgal_count_triangles_in_bbox")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), intent(in) :: bbox_min(3)
      real(c_double), intent(in) :: bbox_max(3)
      integer(c_int) :: count
    end function cgal_count_triangles_in_bbox

    function cgal_collect_triangles_in_bbox(handle, bbox_min, bbox_max, indices, max_count) result(written) &
         bind(C, name="cgal_collect_triangles_in_bbox")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), intent(in) :: bbox_min(3)
      real(c_double), intent(in) :: bbox_max(3)
      integer(c_int), intent(out) :: indices(*)
      integer(c_int), value :: max_count
      integer(c_int) :: written
    end function cgal_collect_triangles_in_bbox
  end interface

contains

  subroutine cgal_initialize_triangle_tree(n_vertices, n_triangles, points, connectivity)
    integer, intent(in) :: n_vertices, n_triangles
    real(c_double), intent(in) :: points(3, n_vertices)
    integer, intent(in) :: connectivity(3, n_triangles)
    real(c_double), allocatable, target :: pts(:, :)
    integer(c_int), allocatable, target :: tris(:, :)
    type(c_ptr) :: handle

    call cgal_finalize_triangle_tree()

    if (n_vertices <= 0 .or. n_triangles <= 0) then
      tree_ready = .false.
      return
    end if

    allocate(pts(3, n_vertices))
    pts = points

    allocate(tris(3, n_triangles))
    tris = int(connectivity, kind=c_int)

    handle = cgal_create_triangle_tree(c_loc(pts(1,1)), int(n_vertices, kind=c_int), &
         c_loc(tris(1,1)), int(n_triangles, kind=c_int))

    if (c_associated(handle)) then
      tree_handle = handle
      tree_ready = .true.
    else
      tree_ready = .false.
    end if

    deallocate(pts, tris)
  end subroutine cgal_initialize_triangle_tree

  subroutine cgal_finalize_triangle_tree()
    if (c_associated(tree_handle)) then
      call cgal_free_triangle_tree(tree_handle)
      tree_handle = c_null_ptr
    end if
    tree_ready = .false.
  end subroutine cgal_finalize_triangle_tree

  logical function cgal_have_triangle_tree()
    cgal_have_triangle_tree = tree_ready
  end function cgal_have_triangle_tree

  subroutine cgal_triangles_in_bbox(bbox_min, bbox_max, indices)
    real(c_double), intent(in) :: bbox_min(3), bbox_max(3)
    integer, allocatable, intent(out) :: indices(:)
    integer(c_int) :: raw_count
    integer :: count
    integer(c_int), allocatable :: buf(:)
    integer :: i

    if (.not. tree_ready) then
      allocate(indices(0))
      return
    end if

    raw_count = cgal_count_triangles_in_bbox(tree_handle, bbox_min, bbox_max)
    count = int(raw_count)
    if (count <= 0) then
      allocate(indices(0))
      return
    end if

    allocate(buf(count))
    raw_count = cgal_collect_triangles_in_bbox(tree_handle, bbox_min, bbox_max, buf, int(count, kind=c_int))
    count = int(raw_count)
    if (count < 0) count = 0
    if (count > size(buf)) count = size(buf)

    if (count <= 0) then
      allocate(indices(0))
      deallocate(buf)
      return
    end if

    allocate(indices(count))
    do i = 1, count
      indices(i) = int(buf(i))
    end do
    deallocate(buf)
  end subroutine cgal_triangles_in_bbox

end module cgal_bindings
