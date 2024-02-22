!....................................................
!..         Subroutine ironing                    ...
!....................................................
subroutine ironingc(f,nx,ny,nz)
implicit none

integer (kind=4),intent(in)   :: nx,ny,nz
complex (kind=8)              :: f(nx,ny,nz)
complex (kind=8), allocatable :: aux(:,:,:)

integer (kind=4)              :: ix,iy,iz
real    (kind=8), parameter   :: a0=8.333333333333333d-2   ! = 1/12
real    (kind=8), parameter   :: a1=0.5d0

allocate(aux(0:(nx+1),0:(ny+1),0:(nz+1)))

aux = 0.d0

ForAll(ix=1:nx,iy=1:ny,iz=1:nz)
  aux(ix,iy,iz) = f(ix,iy,iz)
EndForAll  

aux((nx+1),:,:) = aux(1,:,:)
aux(0,:,:)      = aux(nx,:,:)

aux(:,(ny+1),:) = aux(:,1,:)
aux(:,0,:)      = aux(:,ny,:)

aux(:,:,(nz+1)) = aux(:,:,1)
aux(:,:,0)      = aux(:,:,nz)

!do iz=2,nz-1
!  do iy=2,ny-1
!    do ix=2,nx-1
do iz=1,nz
  do iy=1,ny
    do ix=1,nx
      f(ix,iy,iz) =  a0*( aux(ix-1,iy  ,iz  )+aux(ix+1,iy  ,iz)       &
                         +aux(ix  ,iy-1,iz  )+aux(ix  ,iy+1,iz)       &
                         +aux(ix  ,iy  ,iz-1)+aux(ix  ,iy  ,iz+1) )   &
                  +  a1*aux(ix,iy,iz)
    end do
  end do
end do

deallocate(aux)

return
end
