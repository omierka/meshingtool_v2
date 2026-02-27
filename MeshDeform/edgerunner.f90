  subroutine us_EdgeRunner(f,x,y,z,w,v, ilevel, dcorvg,&
                           kvert,kedge,nel,nvt,net,nat,nProjStep)

  implicit none

  ! Parameters
  real*8, dimension(:) :: f,x,y,z,w,v

  integer, intent(in) :: ilevel

  real*8, dimension(:,:) :: dcorvg

  integer, dimension(:,:), intent(inout) :: kedge, kvert

  integer :: nel,nvt,nat,net

  real*8 :: CHRACTERISTIC_SIZE

  integer :: nProjStep

  ! local variables
  integer NeighE(2,12)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  integer i,j,k,ivt1,ivt2,iProjStep,iaux,iel,iSTL
  real*8 WeightE,P1(3),P2(3),daux2,daux1,PX,PY,PZ,dScale1,dScale2
  real*8 :: dOmega = 0.2d0
  real*8 DIST,dIII
  real*8 :: dCrit1,dCrit2
  real*8 dFactor,dKernel,dPower
  REAL*8, ALLOCATABLE :: myVol(:),DXXX(:),DISTANCE(:)


  ! here we should assign a proper measure to the CHRACTERISTIC_SIZE.
  ! if we use 'box' mesher, lets use the min(delta_x,delta_y,delta_z) from the bounding box
  ! if we use 'annular' mesher, lets use the min(delta_z,0.5*(outer_diam-inner_diam) from the annular bounding box
  CALL SETUP_CHRACTERISTIC_SIZE(CHRACTERISTIC_SIZE)

  ALLOCATE(myVol(nel))
  ALLOCATE(DXXX(nvt))

  if (myid.eq.1) write(*,'(A$)') 'MeshSmoothening:{'
  DO iProjStep=1,nProjStep

   if (myid.eq.1) write(*,'(I0,A1$)') iProjStep," "

   ! here we have to have an immersed call for etimation of the volumes for the elements
   CALL ComputeVolume(myVol)

   f(1:nvt) = 0d0
   w(1:nvt) = 0d0

   DO iel=1,nel
    DO i=1,8
     j = kvert(i,iel)
     f(j) = f(j) + abs(myVol(iel))
     w(j) = w(j) + 1d0
    end DO
   end DO

   DO i=1,nvt
    f(i) = f(i)/w(i)
   end DO

   ! here we have to update the distance values in our vertices against the triangulation and strore that in DISTANCE
   call EstimateDistancesAgainstTriangulation(DISTANCE)

   DO i=1,nvt
     PX = dcorvg(1,i)
     PY = dcorvg(2,i)
     PZ = dcorvg(3,i)

     CALL GetWeight(PX,PY,PZ,dFactor)

     f(i) = dFactor*f(i)
   end DO

   x(1:nvt) = 0d0
   y(1:nvt) = 0d0
   z(1:nvt) = 0d0
   w(1:nvt) = 0d0

   k=1
   DO i=1,nel
    DO j=1,12
     IF (k.eq.kedge(j,i)) THEN
      ivt1 = kvert(NeighE(1,j),i)
      ivt2 = kvert(NeighE(2,j),i)
      P1(:) = dcorvg(:,ivt1)
      P2(:) = dcorvg(:,ivt2)

      daux1 = ABS(f(ivt1))
      daux2 = ABS(f(ivt2))
      WeightE = 1d0

      x(ivt1) = x(ivt1) + WeightE*P2(1)*daux2
      y(ivt1) = y(ivt1) + WeightE*P2(2)*daux2
      z(ivt1) = z(ivt1) + WeightE*P2(3)*daux2
      w(ivt1) = w(ivt1) + WeightE*daux2

      x(ivt2) = x(ivt2) + WeightE*P1(1)*daux1
      y(ivt2) = y(ivt2) + WeightE*P1(2)*daux1
      z(ivt2) = z(ivt2) + WeightE*P1(3)*daux1
      w(ivt2) = w(ivt2) + WeightE*daux1

      k = k + 1
     end IF
    end DO
   end DO

   DO i=1,nvt
     PX = x(i)/w(i)
     PY = y(i)/w(i)
     PZ = z(i)/w(i)
     dcorvg(1,i) = MAX(0d0,(1d0-dOmega))*dcorvg(1,i) + dOmega*PX
     dcorvg(2,i) = MAX(0d0,(1d0-dOmega))*dcorvg(2,i) + dOmega*PY
     dcorvg(3,i) = MAX(0d0,(1d0-dOmega))*dcorvg(3,i) + dOmega*PZ
   end DO

  end DO

  if (myid.eq.1) write(*,'(A1)') "}"

  DEALLOCATE(myVol,DXXX)

   CONTAINS

  subroutine GetWeight(x,y,z,f)
  IMPLICIT NONE
  real*8 x,y,z,f
  real*8 d1
  real*8 f1,dScaleFactor

  dScaleFactor = 100d0/CHRACTERISTIC_SIZE)

  d1 = dScaleFactor*DISTANCE(i)

  CALL KernelFunction(d1,f1)

  f = MIN(f1,25d0)
  f = f**2.3d0

  end IF

  end subroutine GetWeight

  subroutine KernelFunction(d,f)
  real*8 d,daux,f

  if (d.lt.0d0) d = 2.5d0*abs(d)

  IF (d.lt.3.0d0) THEN
   daux = 2d0 + 1.0d0*(3.0d0-d)
  ELSE
   daux = 2d0 - 0.2d0*(d-3.0d0)
  end IF


  f = max(daux,0.8d0)

  end subroutine KernelFunction

  end

