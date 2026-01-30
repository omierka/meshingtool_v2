module tri_tet_intersection
  use, intrinsic :: iso_c_binding, only: c_double
  use cgal_bindings, only: cgal_have_triangle_tree, cgal_triangles_in_bbox
  implicit none
  private
  public :: dp
  public :: tri_tet_intersect, write_intersection_vtu
  public :: tri_hex_intersections, hexa_tet_node_indices
  public :: compute_hex_bounding_sphere, gather_triangles_in_sphere
  public :: compute_hex_subdivision_points
  public :: num_hex_tets, num_hex_sub_vertices

  integer, parameter :: dp = c_double
  integer, parameter :: max_poly_vertices = 32
  integer, parameter :: num_hex_tets = 40
  integer, parameter :: num_hex_corners = 8
  integer, parameter :: num_hex_edges = 12
  integer, parameter :: num_hex_faces = 6
  integer, parameter :: num_hex_sub_vertices = num_hex_corners + num_hex_edges + num_hex_faces + 1
  integer, parameter :: edge_pairs(2, num_hex_edges) = reshape([ &
       1, 2, &
       2, 3, &
       3, 4, &
       4, 1, &
       1, 5, &
       2, 6, &
       3, 7, &
       4, 8, &
       5, 6, &
       6, 7, &
       7, 8, &
       8, 5 ], [2, num_hex_edges])
  integer, parameter :: face_quads(4, num_hex_faces) = reshape([ &
       1, 2, 3, 4, &
       1, 2, 6, 5, &
       2, 3, 7, 6, &
       3, 4, 8, 7, &
       4, 1, 5, 8, &
       5, 6, 7, 8 ], [4, num_hex_faces])
  integer, parameter :: hex_tet_indices(4, num_hex_tets) = reshape([ &
       1,  9, 21, 22, &
       1, 21, 25, 22, &
       1, 21, 12, 25, &
       1, 22, 25, 13, &
      21, 25, 22, 27, &
       2, 10, 21, 23, &
       2, 21, 22, 23, &
       2, 21,  9, 22, &
       2, 23, 22, 14, &
      21, 22, 23, 27, &
       3, 11, 21, 24, &
       3, 21, 23, 24, &
       3, 21, 10, 23, &
       3, 24, 23, 15, &
      21, 23, 24, 27, &
       4, 12, 21, 25, &
       4, 21, 24, 25, &
       4, 21, 11, 24, &
       4, 25, 24, 16, &
      21, 24, 25, 27, &
       5, 17, 26, 22, &
       5, 26, 25, 22, &
       5, 26, 20, 25, &
       5, 22, 25, 13, &
      26, 25, 22, 27, &
       6, 18, 26, 23, &
       6, 26, 22, 23, &
       6, 26, 17, 22, &
       6, 23, 22, 14, &
      26, 22, 23, 27, &
       7, 19, 26, 24, &
       7, 26, 23, 24, &
       7, 26, 18, 23, &
       7, 24, 23, 15, &
      26, 23, 24, 27, &
       8, 20, 26, 25, &
       8, 26, 24, 25, &
       8, 26, 19, 24, &
       8, 25, 24, 16, &
      26, 24, 25, 27 ], [4, num_hex_tets])

contains

  pure function dot3(a,b) result(r)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: r
    r = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end function

  pure function norm2(a) result(r)
    real(dp), intent(in) :: a(3)
    real(dp) :: r
    r = sqrt(dot3(a,a))
  end function

  pure function cross3(a,b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function

  pure function sub3(a,b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = a - b
  end function

  pure function add3(a,b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = a + b
  end function

  pure function mul3(s,a) result(c)
    real(dp), intent(in) :: s
    real(dp), intent(in) :: a(3)
    real(dp) :: c(3)
    c = s*a
  end function

  pure function signed_dist_to_plane(x, n, d) result(sd)
    ! Plane: n·x + d = 0
    real(dp), intent(in) :: x(3), n(3), d
    real(dp) :: sd
    sd = dot3(n,x) + d
  end function

  subroutine make_face_plane(a,b,c, opp, n, d, eps)
    ! Build plane through (a,b,c) with normal oriented so that "inside"
    ! of tetrahedron is sd <= eps (i.e., opposite vertex lies inside).
    real(dp), intent(in) :: a(3), b(3), c(3), opp(3), eps
    real(dp), intent(out) :: n(3), d
    real(dp) :: nn(3), dd, sopp

    nn = cross3(sub3(b,a), sub3(c,a))
    dd = -dot3(nn, a)

    sopp = signed_dist_to_plane(opp, nn, dd)
    ! If opposite vertex is on positive side, flip so that it becomes "inside" (<=0).
    if (sopp > eps) then
      nn = -nn
      dd = -dd
    end if

    n = nn
    d = dd
  end subroutine

  subroutine clip_polygon_by_plane(pin, nin, n, d, eps, pout, nout)
    ! Sutherland–Hodgman polygon clip in 3D against halfspace: sd <= eps
    real(dp), intent(in) :: pin(:, :), n(3), d, eps
    integer, intent(in) :: nin
    real(dp), intent(out) :: pout(:, :)
    integer, intent(out) :: nout

    integer :: i, j
    real(dp) :: s1, s2, t
    real(dp) :: p1(3), p2(3), pi(3)

    nout = 0
    if (nin <= 0) return

    do i = 1, nin
      j = i + 1
      if (j > nin) j = 1

      p1 = pin(:, i)
      p2 = pin(:, j)

      s1 = signed_dist_to_plane(p1, n, d)
      s2 = signed_dist_to_plane(p2, n, d)

      if (s1 <= eps .and. s2 <= eps) then
        ! in->in : keep p2
        nout = nout + 1
        pout(:, nout) = p2

      else if (s1 <= eps .and. s2 > eps) then
        ! in->out : keep intersection
        t = s1 / (s1 - s2)    ! in [0,1]
        pi = add3(p1, mul3(t, sub3(p2,p1)))
        nout = nout + 1
        pout(:, nout) = pi

      else if (s1 > eps .and. s2 <= eps) then
        ! out->in : keep intersection and p2
        t = s1 / (s1 - s2)
        pi = add3(p1, mul3(t, sub3(p2,p1)))
        nout = nout + 1
        pout(:, nout) = pi
        nout = nout + 1
        pout(:, nout) = p2

      else
        ! out->out : keep nothing
      end if
    end do
  end subroutine

  subroutine compact_polygon(p, n, eps)
    ! Remove near-duplicate consecutive points (and closing duplicate)
    real(dp), intent(inout) :: p(:, :)
    integer, intent(inout) :: n
    real(dp), intent(in) :: eps
    integer :: i, k
    real(dp) :: q(3, size(p,2))

    if (n <= 1) return

    k = 0
    do i = 1, n
      if (k == 0) then
        k = 1
        q(:,k) = p(:,i)
      else
        if (norm2(sub3(p(:,i), q(:,k))) > eps) then
          k = k + 1
          q(:,k) = p(:,i)
        end if
      end if
    end do

    ! Remove last if it duplicates first
    if (k >= 2) then
      if (norm2(sub3(q(:,k), q(:,1))) <= eps) k = k - 1
    end if

    p(:,1:k) = q(:,1:k)
    n = k
  end subroutine

  logical pure function is_all_same_point(p, n, eps)
    real(dp), intent(in) :: p(:, :), eps
    integer, intent(in) :: n
    integer :: i
    is_all_same_point = .true.
    if (n <= 1) return
    do i = 2, n
      if (norm2(p(:,i) - p(:,1)) > eps) then
        is_all_same_point = .false.
        return
      end if
    end do
  end function

  logical pure function is_collinear(p, n, eps)
    ! Check if all points lie on a line (within eps)
    real(dp), intent(in) :: p(:, :), eps
    integer, intent(in) :: n
    integer :: i, i2
    real(dp) :: v(3), w(3), c(3)
    is_collinear = .true.
    if (n <= 2) return

    ! Find a direction using first non-identical point
    i2 = 0
    do i = 2, n
      if (norm2(p(:,i) - p(:,1)) > eps) then
        i2 = i
        exit
      end if
    end do
    if (i2 == 0) return  ! all same -> collinear true

    v = p(:,i2) - p(:,1)
    do i = 2, n
      w = p(:,i) - p(:,1)
      c = cross3(v, w)
      if (norm2(c) > eps*max(1.0_dp, norm2(v))) then
        is_collinear = .false.
        return
      end if
    end do
  end function

  subroutine tri_tet_intersect(A,B,C, P,Q,R,S, eps, itype, X, nX)
    ! Compute intersection of triangle ABC with tetrahedron PQRS.
    !
    ! Output:
    !   itype: 0 none, 1 point, 2 segment, 3 polygon
    !   X(:,1:nX): intersection vertices (ordered as clipped polygon)
    !
    real(dp), intent(in) :: A(3),B(3),C(3), P(3),Q(3),R(3),S(3), eps
    integer, intent(out) :: itype, nX
    real(dp), intent(out) :: X(:, :)

    real(dp) :: poly1(3, 32), poly2(3, 32)
    integer :: n1, n2, f
    integer :: i, j, i0, j0
    real(dp) :: npl(3), dpl
    real(dp) :: fa(3), fb(3), fc(3), opp(3)
    real(dp) :: dij, dmax

    ! Start polygon = triangle
    n1 = 3
    poly1(:,1) = A
    poly1(:,2) = B
    poly1(:,3) = C

    ! Clip against the 4 tetrahedron faces.
    do f = 1, 4
      select case (f)
      case (1)
        fa = P; fb = Q; fc = R; opp = S   ! face PQR, opposite S
      case (2)
        fa = P; fb = Q; fc = S; opp = R   ! face PQS, opposite R
      case (3)
        fa = P; fb = R; fc = S; opp = Q   ! face PRS, opposite Q
      case (4)
        fa = Q; fb = R; fc = S; opp = P   ! face QRS, opposite P
      end select

      call make_face_plane(fa,fb,fc, opp, npl, dpl, eps)
      call clip_polygon_by_plane(poly1, n1, npl, dpl, eps, poly2, n2)

      n1 = n2
      if (n1 <= 0) exit
      poly1(:,1:n1) = poly2(:,1:n1)
      call compact_polygon(poly1, n1, eps)
      if (n1 <= 0) exit
    end do

    ! Copy to output
    nX = n1
    if (nX > 0) X(:,1:nX) = poly1(:,1:nX)

    ! Classify
    if (nX == 0) then
      itype = 0
      return
    end if

    call compact_polygon(X, nX, eps)

    if (nX == 0) then
      itype = 0
      return
    end if

    if (is_all_same_point(X, nX, eps)) then
      itype = 1
      nX = 1
      return
    end if

    if (nX == 2) then
      if (norm2(X(:,2)-X(:,1)) <= eps) then
        itype = 1
        nX = 1
      else
        itype = 2
      end if
      return
    end if

    if (is_collinear(X, nX, eps)) then
      ! Reduce to segment endpoints (simple approach: pick farthest pair)
      dmax = -1.0_dp; i0 = 1; j0 = 2
      do i = 1, nX
        do j = i+1, nX
          dij = norm2(X(:,j)-X(:,i))
          if (dij > dmax) then
            dmax = dij
            i0 = i; j0 = j
          end if
        end do
      end do
      if (dmax <= eps) then
        itype = 1
        X(:,1) = X(:,i0)
        nX = 1
      else
        itype = 2
        X(:,1) = X(:,i0)
        X(:,2) = X(:,j0)
        nX = 2
      end if
      return
    end if

    itype = 3
    ! nX >= 3: polygon vertices already ordered by clipping
  end subroutine

  subroutine hexa_tet_node_indices(hexIndices, tetNodes)
    integer, intent(in) :: hexIndices(num_hex_sub_vertices)
    integer, intent(out) :: tetNodes(4, num_hex_tets)
    integer :: t
    do t = 1, num_hex_tets
      tetNodes(:, t) = hexIndices(hex_tet_indices(:, t))
    end do
  end subroutine hexa_tet_node_indices

  pure subroutine compute_hex_subdivision_points(hexPts, subdivPts)
    real(dp), intent(in) :: hexPts(3, num_hex_corners)
    real(dp), intent(out) :: subdivPts(3, num_hex_sub_vertices)
    integer :: e, f, idx
    real(dp) :: inv_four

    subdivPts(:, 1:num_hex_corners) = hexPts

    do e = 1, num_hex_edges
      idx = num_hex_corners + e
      subdivPts(:, idx) = 0.5_dp * (hexPts(:, edge_pairs(1,e)) + hexPts(:, edge_pairs(2,e)))
    end do

    inv_four = 0.25_dp
    do f = 1, num_hex_faces
      idx = num_hex_corners + num_hex_edges + f
      subdivPts(:, idx) = inv_four * ( hexPts(:, face_quads(1,f)) + hexPts(:, face_quads(2,f)) + &
                                       hexPts(:, face_quads(3,f)) + hexPts(:, face_quads(4,f)) )
    end do

    subdivPts(:, num_hex_sub_vertices) = sum(hexPts(:,1:num_hex_corners), dim=2) / real(num_hex_corners, dp)
  end subroutine compute_hex_subdivision_points

  pure subroutine hexa_tet_vertices(hexPts, tetVerts)
    real(dp), intent(in) :: hexPts(3,8)
    real(dp), intent(out) :: tetVerts(3,4,num_hex_tets)
    real(dp) :: subdivPts(3, num_hex_sub_vertices)
    integer :: t, k

    call compute_hex_subdivision_points(hexPts, subdivPts)

    do t = 1, num_hex_tets
      do k = 1, 4
        tetVerts(:,k,t) = subdivPts(:, hex_tet_indices(k,t))
      end do
    end do
  end subroutine hexa_tet_vertices

  subroutine compute_hex_bounding_sphere(hexPts, center, radius)
    real(dp), intent(in) :: hexPts(3,8)
    real(dp), intent(out) :: center(3)
    real(dp), intent(out) :: radius

    call minimal_enclosing_sphere(hexPts, 8, center, radius)
  end subroutine compute_hex_bounding_sphere

  subroutine gather_triangles_in_sphere(center, radius, hexPts, tri_points, tri_connectivity, tri_count, indices)
    real(dp), intent(in) :: center(3), radius
    real(dp), intent(in) :: hexPts(3,8)
    real(dp), intent(in) :: tri_points(:, :)
    integer, intent(in) :: tri_connectivity(:, :)
    integer, intent(in) :: tri_count
    integer, allocatable, intent(out) :: indices(:)
    real(dp) :: bbox_min(3), bbox_max(3)
    integer, allocatable :: cgal_candidates(:)
    integer, allocatable :: selected(:)
    integer :: i, tri_id, count
    real(dp) :: triA(3), triB(3), triC(3)
    logical :: include_triangle

    call compute_hex_bbox(hexPts, bbox_min, bbox_max)
    if (.not. cgal_have_triangle_tree()) then
      allocate(indices(0))
      return
    end if

    call cgal_triangles_in_bbox(bbox_min, bbox_max, cgal_candidates)
    if (size(cgal_candidates) > 0) then
      allocate(selected(size(cgal_candidates)))
      count = 0
      do i = 1, size(cgal_candidates)
        tri_id = cgal_candidates(i)
        if (tri_id < 1 .or. tri_id > tri_count) cycle
        triA = tri_points(:, tri_connectivity(1, tri_id))
        triB = tri_points(:, tri_connectivity(2, tri_id))
        triC = tri_points(:, tri_connectivity(3, tri_id))
        include_triangle = triangle_overlaps_bbox(bbox_min, bbox_max, triA, triB, triC) .or. &
             triangle_vertex_inside_bbox(bbox_min, bbox_max, triA, triB, triC)
        if (include_triangle) then
          count = count + 1
          selected(count) = tri_id
        end if
      end do
      if (count > 0) then
        allocate(indices(count))
        indices = selected(1:count)
      else
        allocate(indices(0))
      end if
      deallocate(selected)
    else
      allocate(indices(0))
    end if
    if (allocated(cgal_candidates)) deallocate(cgal_candidates)

  end subroutine gather_triangles_in_sphere

  subroutine gather_triangles_in_sphere_fallback(center, radius, hexPts, tri_points, tri_connectivity, tri_count, indices)
    real(dp), intent(in) :: center(3), radius
    real(dp), intent(in) :: hexPts(3,8)
    real(dp), intent(in) :: tri_points(:, :)
    integer, intent(in) :: tri_connectivity(:, :)
    integer, intent(in) :: tri_count
    integer, allocatable, intent(out) :: indices(:)

    integer :: i, count
    logical :: inside

    count = 0
    do i = 1, tri_count
      inside = triangle_intersects_sphere(center, radius, &
           tri_points(:, tri_connectivity(1,i)), &
           tri_points(:, tri_connectivity(2,i)), &
           tri_points(:, tri_connectivity(3,i))) .or. &
           any_triangle_vertex_inside_hex(hexPts, &
           tri_points(:, tri_connectivity(1,i)), &
           tri_points(:, tri_connectivity(2,i)), &
           tri_points(:, tri_connectivity(3,i)))
      if (inside) count = count + 1
    end do

    if (count <= 0) then
      allocate(indices(0))
      return
    end if

    allocate(indices(count))
    count = 0
    do i = 1, tri_count
      inside = triangle_intersects_sphere(center, radius, &
           tri_points(:, tri_connectivity(1,i)), &
           tri_points(:, tri_connectivity(2,i)), &
           tri_points(:, tri_connectivity(3,i))) .or. &
           any_triangle_vertex_inside_hex(hexPts, &
           tri_points(:, tri_connectivity(1,i)), &
           tri_points(:, tri_connectivity(2,i)), &
           tri_points(:, tri_connectivity(3,i)))
      if (inside) then
        count = count + 1
        indices(count) = i
      end if
    end do
  end subroutine gather_triangles_in_sphere_fallback

  pure subroutine compute_hex_bbox(hexPts, bbox_min, bbox_max)
    real(dp), intent(in) :: hexPts(3,8)
    real(dp), intent(out) :: bbox_min(3), bbox_max(3)
    integer :: k
    do k = 1, 3
      bbox_min(k) = minval(hexPts(k, :))
      bbox_max(k) = maxval(hexPts(k, :))
    end do
  end subroutine compute_hex_bbox

  pure logical function boxes_overlap(minA, maxA, minB, maxB)
    real(dp), intent(in) :: minA(3), maxA(3), minB(3), maxB(3)
    integer :: i
    boxes_overlap = .true.
    do i = 1, 3
      if (maxA(i) < minB(i)) then
        boxes_overlap = .false.
        return
      end if
      if (minA(i) > maxB(i)) then
        boxes_overlap = .false.
        return
      end if
    end do
  end function boxes_overlap

  pure logical function point_inside_bbox(bbox_min, bbox_max, p)
    real(dp), intent(in) :: bbox_min(3), bbox_max(3), p(3)
    real(dp), parameter :: tol = 1.0d-12
    integer :: i
    point_inside_bbox = .true.
    do i = 1, 3
      if (p(i) < bbox_min(i) - tol) then
        point_inside_bbox = .false.
        return
      end if
      if (p(i) > bbox_max(i) + tol) then
        point_inside_bbox = .false.
        return
      end if
    end do
  end function point_inside_bbox

  pure logical function triangle_vertex_inside_bbox(bbox_min, bbox_max, a, b, c)
    real(dp), intent(in) :: bbox_min(3), bbox_max(3)
    real(dp), intent(in) :: a(3), b(3), c(3)
    triangle_vertex_inside_bbox = point_inside_bbox(bbox_min, bbox_max, a) .or. &
         point_inside_bbox(bbox_min, bbox_max, b) .or. &
         point_inside_bbox(bbox_min, bbox_max, c)
  end function triangle_vertex_inside_bbox

  pure logical function triangle_overlaps_bbox(bbox_min, bbox_max, a, b, c)
    real(dp), intent(in) :: bbox_min(3), bbox_max(3)
    real(dp), intent(in) :: a(3), b(3), c(3)
    real(dp) :: tri_min(3), tri_max(3)
    integer :: k
    do k = 1, 3
      tri_min(k) = minval([a(k), b(k), c(k)])
      tri_max(k) = maxval([a(k), b(k), c(k)])
    end do
    triangle_overlaps_bbox = boxes_overlap(bbox_min, bbox_max, tri_min, tri_max)
  end function triangle_overlaps_bbox


  pure logical function any_triangle_vertex_inside_hex(hexPts, p1, p2, p3)
    real(dp), intent(in) :: hexPts(3,8)
    real(dp), intent(in) :: p1(3), p2(3), p3(3)
    real(dp) :: tetVerts(3,4,num_hex_tets)
    integer :: t
    call hexa_tet_vertices(hexPts, tetVerts)
    any_triangle_vertex_inside_hex = .false.
    do t = 1, num_hex_tets
      if (point_in_tetra(p1, tetVerts(:,1,t), tetVerts(:,2,t), tetVerts(:,3,t), tetVerts(:,4,t))) then
        any_triangle_vertex_inside_hex = .true.
        return
      end if
      if (point_in_tetra(p2, tetVerts(:,1,t), tetVerts(:,2,t), tetVerts(:,3,t), tetVerts(:,4,t))) then
        any_triangle_vertex_inside_hex = .true.
        return
      end if
      if (point_in_tetra(p3, tetVerts(:,1,t), tetVerts(:,2,t), tetVerts(:,3,t), tetVerts(:,4,t))) then
        any_triangle_vertex_inside_hex = .true.
        return
      end if
    end do
  end function any_triangle_vertex_inside_hex

  pure logical function point_in_tetra(p, a, b, c, d)
    real(dp), intent(in) :: p(3), a(3), b(3), c(3), d(3)
    real(dp) :: mat(3,3), rhs(3), bary(3)
    logical :: success
    real(dp) :: w0, eps
    eps = 1.0e-10_dp
    mat(:,1) = b - a
    mat(:,2) = c - a
    mat(:,3) = d - a
    rhs = p - a
    call solve_linear3(mat, rhs, bary, success)
    if (.not. success) then
      point_in_tetra = .false.
      return
    end if
    w0 = 1.0_dp - bary(1) - bary(2) - bary(3)
    point_in_tetra = (bary(1) >= -eps .and. bary(2) >= -eps .and. bary(3) >= -eps .and. w0 >= -eps .and. &
                      bary(1) <= 1.0_dp + eps .and. bary(2) <= 1.0_dp + eps .and. bary(3) <= 1.0_dp + eps .and. &
                      w0 <= 1.0_dp + eps)
  end function point_in_tetra

  pure logical function triangle_intersects_sphere(center, radius, p1, p2, p3)
    real(dp), intent(in) :: center(3), radius
    real(dp), intent(in) :: p1(3), p2(3), p3(3)
    real(dp) :: dist2, rsq, eps
    eps = 1.0e-12_dp
    rsq = radius*radius
    dist2 = point_triangle_distance2(center, p1, p2, p3)
    triangle_intersects_sphere = (dist2 <= rsq + eps)
  end function triangle_intersects_sphere

  pure real(dp) function point_triangle_distance2(p, a, b, c)
    real(dp), intent(in) :: p(3), a(3), b(3), c(3)
    real(dp) :: ab(3), ac(3), ap(3), bp(3), cp(3), bc(3)
    real(dp) :: d1, d2, d3, d4, d5, d6, vc, vb, va, v, w
    real(dp) :: d00, d01, d02, d11, d12, denom

    ab = b - a
    ac = c - a
    ap = p - a
    d1 = dot3(ab, ap)
    d2 = dot3(ac, ap)
    if (d1 <= 0.0_dp .and. d2 <= 0.0_dp) then
      point_triangle_distance2 = dot3(ap, ap)
      return
    end if

    bp = p - b
    d3 = dot3(ab, bp)
    d4 = dot3(ac, bp)
    if (d3 >= 0.0_dp .and. d4 <= d3) then
      point_triangle_distance2 = dot3(bp, bp)
      return
    end if

    vc = d1*d4 - d3*d2
    if (vc <= 0.0_dp .and. d1 >= 0.0_dp .and. d3 <= 0.0_dp) then
      v = d1 / (d1 - d3)
      point_triangle_distance2 = dot3(p - (a + v*ab), p - (a + v*ab))
      return
    end if

    cp = p - c
    d5 = dot3(ab, cp)
    d6 = dot3(ac, cp)
    if (d6 >= 0.0_dp .and. d5 <= d6) then
      point_triangle_distance2 = dot3(cp, cp)
      return
    end if

    vb = d5*d2 - d1*d6
    if (vb <= 0.0_dp .and. d2 >= 0.0_dp .and. d6 <= 0.0_dp) then
      w = d2 / (d2 - d6)
      point_triangle_distance2 = dot3(p - (a + w*ac), p - (a + w*ac))
      return
    end if

    bc = c - b
    va = d3*d6 - d5*d4
    if (va <= 0.0_dp .and. (d4 - d3) >= 0.0_dp .and. (d5 - d6) >= 0.0_dp) then
      w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
      point_triangle_distance2 = dot3(p - (b + w*bc), p - (b + w*bc))
      return
    end if

    d00 = dot3(ac, ac)
    d01 = dot3(ac, ab)
    d02 = dot3(ac, ap)
    d11 = dot3(ab, ab)
    d12 = dot3(ab, ap)
    denom = d00*d11 - d01*d01
    if (denom <= 1.0e-16_dp) then
      point_triangle_distance2 = dot3(ap, ap)
      return
    end if
    denom = 1.0_dp / denom
    v = (d11*d02 - d01*d12) * denom
    w = (d00*d12 - d01*d02) * denom
    point_triangle_distance2 = dot3(p - (a + v*ac + w*ab), p - (a + v*ac + w*ab))
  end function point_triangle_distance2

  subroutine minimal_enclosing_sphere(points, npts, center, radius)
    real(dp), intent(in) :: points(:, :)
    integer, intent(in) :: npts
    real(dp), intent(out) :: center(3)
    real(dp), intent(out) :: radius

    real(dp) :: best_center(3), temp_center(3)
    real(dp) :: best_radius, temp_radius
    real(dp), parameter :: eps = 1.0e-12_dp
    integer :: i, j, k, l
    logical :: success, contains, found

    best_radius = huge(1.0_dp)
    found = .false.

    if (npts <= 0) then
      center = 0.0_dp
      radius = 0.0_dp
      return
    end if

    ! Single points
    do i = 1, npts
      temp_center = points(:, i)
      temp_radius = 0.0_dp
      contains = all_points_within(points, npts, temp_center, temp_radius, eps)
      if (contains) then
        best_center = temp_center
        best_radius = temp_radius
        found = .true.
        exit
      end if
    end do

    ! Pairs
    do i = 1, npts
      do j = i+1, npts
        call sphere_from_pair(points(:,i), points(:,j), temp_center, temp_radius)
        contains = all_points_within(points, npts, temp_center, temp_radius, eps)
        if (contains) then
          if (.not. found .or. temp_radius < best_radius) then
            best_center = temp_center
            best_radius = temp_radius
            found = .true.
          end if
        end if
      end do
    end do

    ! Triples
    do i = 1, npts
      do j = i+1, npts
        do k = j+1, npts
          call sphere_from_triplet(points(:,i), points(:,j), points(:,k), temp_center, temp_radius, success)
          if (.not. success) cycle
          contains = all_points_within(points, npts, temp_center, temp_radius, eps)
          if (contains) then
            if (.not. found .or. temp_radius < best_radius) then
              best_center = temp_center
              best_radius = temp_radius
              found = .true.
            end if
          end if
        end do
      end do
    end do

    ! Quadruples
    do i = 1, npts
      do j = i+1, npts
        do k = j+1, npts
          do l = k+1, npts
            call sphere_from_quad(points(:,i), points(:,j), points(:,k), points(:,l), temp_center, temp_radius, success)
            if (.not. success) cycle
            contains = all_points_within(points, npts, temp_center, temp_radius, eps)
            if (contains) then
              if (.not. found .or. temp_radius < best_radius) then
                best_center = temp_center
                best_radius = temp_radius
                found = .true.
              end if
            end if
          end do
        end do
      end do
    end do

    if (found) then
      center = best_center
      radius = best_radius
    else
      call bounding_box_sphere(points, npts, center, radius)
    end if
  end subroutine minimal_enclosing_sphere

  logical function all_points_within(points, npts, center, radius, eps)
    real(dp), intent(in) :: points(:, :)
    integer, intent(in) :: npts
    real(dp), intent(in) :: center(3), radius, eps
    integer :: i
    all_points_within = .true.
    do i = 1, npts
      if (norm2(points(:,i) - center) > radius + eps) then
        all_points_within = .false.
        return
      end if
    end do
  end function all_points_within

  subroutine sphere_from_pair(p1, p2, center, radius)
    real(dp), intent(in) :: p1(3), p2(3)
    real(dp), intent(out) :: center(3), radius
    center = 0.5_dp * (p1 + p2)
    radius = 0.5_dp * norm2(p1 - p2)
  end subroutine sphere_from_pair

  subroutine sphere_from_triplet(pa, pb, pc, center, radius, success)
    real(dp), intent(in) :: pa(3), pb(3), pc(3)
    real(dp), intent(out) :: center(3), radius
    logical, intent(out) :: success
    real(dp) :: matA(3,3), rhsvec(3), normal(3)
    real(dp) :: xvec(3)
    success = .false.

    normal = cross3(pb - pa, pc - pa)
    if (norm2(normal) < 1.0e-12_dp) return

    matA(1,:) = pb - pa
    matA(2,:) = pc - pa
    matA(3,:) = normal
    rhsvec(1) = 0.5_dp*(dot3(pb,pb) - dot3(pa,pa))
    rhsvec(2) = 0.5_dp*(dot3(pc,pc) - dot3(pa,pa))
    rhsvec(3) = dot3(normal, pa)

    call solve_linear3(matA, rhsvec, xvec, success)
    if (.not. success) return

    center = xvec
    radius = norm2(center - pa)
  end subroutine sphere_from_triplet

  subroutine sphere_from_quad(pa, pb, pc, pd, center, radius, success)
    real(dp), intent(in) :: pa(3), pb(3), pc(3), pd(3)
    real(dp), intent(out) :: center(3), radius
    logical, intent(out) :: success
    real(dp) :: matA(3,3), rhsvec(3), xvec(3)
    success = .false.

    matA(1,:) = pa - pd
    matA(2,:) = pb - pd
    matA(3,:) = pc - pd
    rhsvec(1) = 0.5_dp*(dot3(pa,pa) - dot3(pd,pd))
    rhsvec(2) = 0.5_dp*(dot3(pb,pb) - dot3(pd,pd))
    rhsvec(3) = 0.5_dp*(dot3(pc,pc) - dot3(pd,pd))

    call solve_linear3(matA, rhsvec, xvec, success)
    if (.not. success) return

    center = xvec
    radius = norm2(center - pa)
  end subroutine sphere_from_quad

  subroutine bounding_box_sphere(points, npts, center, radius)
    real(dp), intent(in) :: points(:, :)
    integer, intent(in) :: npts
    real(dp), intent(out) :: center(3), radius
    real(dp) :: minv(3), maxv(3)
    integer :: i
    minv = points(:,1)
    maxv = points(:,1)
    do i = 2, npts
      minv = min(minv, points(:,i))
      maxv = max(maxv, points(:,i))
    end do
    center = 0.5_dp*(minv + maxv)
    radius = norm2(maxv - center)
  end subroutine bounding_box_sphere

  pure subroutine solve_linear3(Ain, bin, x, success)
    real(dp), intent(in) :: Ain(3,3)
    real(dp), intent(in) :: bin(3)
    real(dp), intent(out) :: x(3)
    logical, intent(out) :: success
    real(dp) :: A(3,3), bvec(3), factor
    integer :: i, j, pivot
    real(dp), parameter :: eps = 1.0e-12_dp

    A = Ain
    bvec = bin
    success = .false.

    do i = 1, 3
      pivot = i
      do j = i+1, 3
        if (abs(A(j,i)) > abs(A(pivot,i))) pivot = j
      end do
      if (abs(A(pivot,i)) < eps) return
      if (pivot /= i) then
        call swap_rows(A, bvec, i, pivot)
      end if
      factor = A(i,i)
      A(i,:) = A(i,:) / factor
      bvec(i) = bvec(i) / factor
      do j = 1, 3
        if (j == i) cycle
        factor = A(j,i)
        A(j,:) = A(j,:) - factor*A(i,:)
        bvec(j) = bvec(j) - factor*bvec(i)
      end do
    end do

    x = bvec
    success = .true.
  end subroutine solve_linear3

  pure subroutine swap_rows(A, bvec, i, j)
    real(dp), intent(inout) :: A(3,3)
    real(dp), intent(inout) :: bvec(3)
    integer, intent(in) :: i, j
    real(dp) :: temp
    real(dp) :: rowtemp(3)

    rowtemp = A(i,:)
    A(i,:) = A(j,:)
    A(j,:) = rowtemp

    temp = bvec(i)
    bvec(i) = bvec(j)
    bvec(j) = temp
  end subroutine swap_rows

  subroutine tri_hex_intersections(A,B,C, hexPts, eps, itypes, nXs, Xs)
    real(dp), intent(in) :: A(3), B(3), C(3)
    real(dp), intent(in) :: hexPts(3,8)
    real(dp), intent(in) :: eps
    integer, intent(out) :: itypes(num_hex_tets)
    integer, intent(out) :: nXs(num_hex_tets)
    real(dp), intent(out) :: Xs(3, max_poly_vertices, num_hex_tets)

    real(dp) :: tetVerts(3,4,num_hex_tets)
    integer :: t

    call hexa_tet_vertices(hexPts, tetVerts)

    do t = 1, num_hex_tets
      call tri_tet_intersect(A,B,C, tetVerts(:,1,t), tetVerts(:,2,t), tetVerts(:,3,t), &
                             tetVerts(:,4,t), eps, itypes(t), Xs(:,:,t), nXs(t))
    end do
  end subroutine tri_hex_intersections

  subroutine write_intersection_vtu(filename, A,B,C, P,Q,R,S, X, nX, itype)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: A(3), B(3), C(3)
    real(dp), intent(in) :: P(3), Q(3), R(3), S(3)
    real(dp), intent(in) :: X(:, :)
    integer, intent(in) :: nX, itype

    integer :: unit, ios
    integer :: nPoints, nCells, baseInt
    integer :: nconn_int, cellTypeInt
    integer :: total_conn
    integer, allocatable :: conn(:), offsets(:), cellTypes(:), geom_ids(:)
    integer :: i

    nPoints = 7
    if (nX > 0) nPoints = nPoints + nX

    nCells = 2
    nconn_int = 0
    cellTypeInt = 0
    if (itype > 0 .and. nX > 0) then
      nCells = 3
      select case (itype)
      case (1)
        nconn_int = 1
        cellTypeInt = 1          ! VTK_VERTEX
      case (2)
        nconn_int = 2
        cellTypeInt = 3          ! VTK_LINE
      case default
        nconn_int = nX
        cellTypeInt = 7          ! VTK_POLYGON
      end select
    end if

    baseInt = 7

    allocate(offsets(nCells), cellTypes(nCells), geom_ids(nCells))
    offsets(1) = 4
    cellTypes(1) = 10            ! VTK_TETRA
    geom_ids(1) = 0
    offsets(2) = 7
    cellTypes(2) = 5             ! VTK_TRIANGLE
    geom_ids(2) = 1
    if (nCells == 3) then
      offsets(3) = offsets(2) + nconn_int
      cellTypes(3) = cellTypeInt
      geom_ids(3) = 2
    end if
    total_conn = offsets(nCells)

    allocate(conn(total_conn))
    conn(1:4) = [0,1,2,3]
    conn(5:7) = [4,5,6]
    if (nCells == 3) then
      do i = 1, nconn_int
        conn(7+i) = baseInt + i - 1
      end do
    end if

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Failed to open file ', trim(filename)
      deallocate(conn, offsets, cellTypes, geom_ids)
      return
    end if

    write(unit,'(A)') '<?xml version="1.0"?>'
    write(unit,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(unit,'(A)') '  <UnstructuredGrid>'
    write(unit,'(A," NumberOfPoints=""",I0,""" NumberOfCells=""",I0,""">")') '    <Piece', nPoints, nCells
    write(unit,'(A)') '      <PointData/>'
    write(unit,'(A)') '      <CellData Scalars="geometry_id">'
    write(unit,'(A)') '        <DataArray type="Int32" Name="geometry_id" format="ascii">'
    write(unit,'(A,*(I0,1X))') '          ', geom_ids
    write(unit,'(A)') '        </DataArray>'
    write(unit,'(A)') '      </CellData>'
    write(unit,'(A)') '      <Points>'
    write(unit,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    call write_point_row(P)
    call write_point_row(Q)
    call write_point_row(R)
    call write_point_row(S)
    call write_point_row(A)
    call write_point_row(B)
    call write_point_row(C)
    if (nX > 0) then
      do i = 1, nX
        call write_point_row(X(:,i))
      end do
    end if
    write(unit,'(A)') '        </DataArray>'
    write(unit,'(A)') '      </Points>'
    write(unit,'(A)') '      <Cells>'
    write(unit,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
    write(unit,'(A,*(I0,1X))') '          ', conn
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

    deallocate(conn, offsets, cellTypes, geom_ids)

  contains

    subroutine write_point_row(pt)
      real(dp), intent(in) :: pt(3)
      write(unit,'(A,3(1X,ES23.15))') '          ', pt(1), pt(2), pt(3)
    end subroutine write_point_row

  end subroutine write_intersection_vtu

end module tri_tet_intersection
