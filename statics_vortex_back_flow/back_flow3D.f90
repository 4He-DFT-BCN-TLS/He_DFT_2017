
subroutine energy_v
Use fftmodule
Use backflow
use alphasterm
use deriva
use energies
use field
use lenard4
use grid
use gridk
use he4
use rho
use util1
use work1

implicit none
real (kind=8)              ::  aux1,aux2,aux3,aux4,aux5,fac,amass ! Auxiliar variables
complex (kind=8)           ::  ci,auxc
integer (kind=4) :: ix,iy,iz

! pot_bf_aux d'ha d'utilitzar si es vol fer un cutoff suau

ci=(0.d0,1.d0)
fac=0.0763822d0
amass=0.000481402d0

! It must computed before:  fden

! Now compute:  integral (dr integral( dr' U(r-r')*rho(r') rho(r) v(r)**2 )

wk1 = fden*ulj
call fftbk_1()  ! wk1 |----> sto1

! This is v(r)**2

Sto5=vel(1,:,:,:)**2+vel(2,:,:,:)**2+vel(3,:,:,:)**2

sto4 = sto1*sto5*den

!pot_bf_aux = sto4

aux1  =sum(sto4) *dxyz

!--------------------------------------------------------------------
! Now compute:  integral (dr integral( dr' U(r-r')*rho(r') rho(r) v(r')**2 )

sto1 = sto5*den

call fftfw_1()  ! sto1 |----> wk1
wk1 = wk1*ulj
call fftbk_1()  !  wk1 |----> sto1
sto4 = sto1*den

!pot_bf_aux = pot_bf_aux + sto4

aux2  =sum(sto4) *dxyz

!--------------------------------------------------------------------
! Now compute:  -2* integral (dr integral( dr' U(r-r')*rho(r') rho(r) v_x(r)*v_x(r') )
!               -2* integral (dr integral( dr' U(r-r')*rho(r') rho(r) v_y(r)*v_y(r') )
!               -2* integral (dr integral( dr' U(r-r')*rho(r') rho(r) v_z(r)*v_z(r') )

sto1=den*vel(1,:,:,:)
sto2=den*vel(2,:,:,:)
sto3=den*vel(3,:,:,:)

call fftfw_123()    ! sto1, 2, 3 |----> wk1, 2, 3
wk1 = wk1*ulj
wk2 = wk2*ulj
wk3 = wk3*ulj
call fftbk_xyz()    !  wk1, 2, 3 |----> intxalf, intyalf, intyalf

sto4=den*(intxalf*vel(1,:,:,:) + intyalf*vel(2,:,:,:) + intzalf*vel(3,:,:,:))

!pot_bf_aux = pot_bf_aux -2.*sto4

aux3  =-2.*sum(sto4) *dxyz


! This is E_bf

!Write(6,'(/,"Back flow energy terms:",1p,5E15.6,/)')aux1,aux2,aux3



e_bf=-0.25d0*amass*(aux1+aux2+aux3)

!pot_bf_aux = -0.25d0*amass*pot_bf_aux


return

end


subroutine poten_v
Use fftmodule
Use Para_DerivnD
Use backflow
use alphasterm
use deriva
use energies
use field
use lenard4
use grid
use gridk
use he4
use rho
use util1
use work1

implicit none
real (kind=8)              ::  aux1,aux2,aux3,aux4,aux5,fac,amass ! Auxiliar variables
complex (kind=8)           ::  ci,auxc
integer (kind=4) :: ix,iy,iz,is


ci=(0.d0,1.d0)
fac=0.0763822d0
amass=0.000481402d0

call velocity()

! It must computed before:  fden, dxden, dyden, dzden


! Now compute:  integral( dr' U(r-r')*rho(r') )

!Write(6,'("From poten_v...(1)")')
!Call flush(6)

wk1 = fden*ulj  
call fftbk_4() !  wk1 |----> sto4
!Write(6,'("From poten_v...(2)")')
!Call flush(6)


! This is v(r)**2

sto5=vel(1,:,:,:)**2+vel(2,:,:,:)**2+vel(3,:,:,:)**2

sto2 = sto4*sto5

!--------------------------------------------------------------------
! Now compute:  integral( dr' U(r-r')*rho(r') v(r')**2 )

sto1 = sto5*den
!Write(6,'("From poten_v...(Previ-3)")')
!Call flush(6)
call fftfw_1() ! sto1 |----> wk1
!Write(6,'("From poten_v...(3)")')
!Call flush(6)

wk1 = wk1*ulj  
call fftbk_1() !  wk1 |---> sto1  
!Write(6,'("From poten_v...(4)")')
!Call flush(6)

sto2 = sto2 + sto1

!--------------------------------------------------------------------
! Now compute:  -2* integral( dr' U(r-r')*rho(r')  v_x(r)*v_x(r') )
!               -2* integral( dr' U(r-r')*rho(r')  v_y(r)*v_y(r') )
!               -2* integral( dr' U(r-r')*rho(r')  v_z(r)*v_z(r') )

sto1=den*vel(1,:,:,:)
sto2=den*vel(2,:,:,:)
sto3=den*vel(3,:,:,:)

call fftfw_123()    ! sto1, 2, 3 |----> wk1, wk2, wk3
!Write(6,'("From poten_v...(5)")')
!Call flush(6)

wk1 = wk1*ulj
wk2 = wk2*ulj
wk3 = wk3*ulj
call fftbk_xyz()    !  wk1, 2, 3 |----> intxalf, intyalf, intyalf
!Write(6,'("From poten_v...(6)")')
!Call flush(6)


!--------------------------------------------------------------------
! Now compute:  -2* integral( dr' U(r-r')*rho(r')  v_x(r)*v_x(r') ) +(..y..) +(..z..)
!--------------------------------------------------------------------

sto2 = sto2 -2.*(intxalf*vel(1,:,:,:) + intyalf*vel(2,:,:,:) + intzalf*vel(3,:,:,:))
!Write(6,'("From poten_v...(7)")')
!Call flush(6)

! hydrodynamic contribution to the potential
pot_bf =  - 0.5*amass*sto2

! non local contribution

!--------------------------
! Now compute:  (grad(ro).v integral( dr' U(r-r')*rho(r') )
!--------------------------


!sto6=dxden*vel(1,:,:,:)
!sto6=sto6+dyden*vel(2,:,:,:)
!sto6=sto6+dzden*vel(3,:,:,:)
!sto5=sto6*sto4/den

!--------------------------
! Now compute:  (grad(ro).v integral( dr' U(r-r')*rho(r') )
!--------------------------
sto5=(dxden*vel(1,:,:,:)+dyden*vel(2,:,:,:)+dzden*vel(3,:,:,:))*sto4/den    &
!--------------------------
! Now compute:  -grad(ro)/ro integral( dr' U(r-r')*rho(r') v(r') )
!--------------------------
    -(dxden*intxalf+dyden*intyalf+dzden*intzalf)/den 
!Write(6,'("From poten_v...(8)")')
!Call flush(6)

!sto5=sto5-(dxden*intxalf+dyden*intyalf+dzden*intzalf)/den

!--------------------------
! Now compute:  div(v) integral( dr' U(r-r')*rho(r') )
!--------------------------
!Write(6,'("From poten_v...(Previ-9)")')
!Call flush(6)

sto2 = vel(1,:,:,:)
Call derivnD(1,nn,hx,1,sto2,sto1,Icon)
!Call derivnD(1,nn,hx,1,vel(1,:,:,:),sto1,Icon)
sto5=sto5+sto1*sto4

sto2 = vel(2,:,:,:)
Call derivnD(1,nn,hy,2,sto2,sto1,Icon)
!Call derivnD(1,nn,hy,2,vel(2,:,:,:),sto1,Icon)
sto5=sto5+sto1*sto4

sto2 = vel(3,:,:,:)
Call derivnD(1,nn,hz,3,sto2,sto1,Icon)
!Call derivnD(1,nn,hz,3,vel(3,:,:,:),sto1,Icon)
sto5=sto5+sto1*sto4

!Write(6,'("From poten_v...(9)")')
!Call flush(6)

!--------------------------------------------------------------------
! Now compute:  - div (ntegral( dr' U(r-r')*rho(r') v(r') )
!--------------------------------------------------------------------

Call derivnD(1,nn,hx,1,intxalf,sto1,Icon)
sto5=sto5-sto1

Call derivnD(1,nn,hy,2,intyalf,sto1,Icon)
sto5=sto5-sto1

Call derivnD(1,nn,hz,3,intzalf,sto1,Icon)
sto5=sto5-sto1

!Write(6,'("From poten_v...(10)")')
!Call flush(6)

pot_bf =  pot_bf + 0.5*fac*ci*sto5

return

end

subroutine velocity

Use Para_DerivnD
Use backflow
use deriva ! (icon,npd,dxden,dyden,dzden,pderx,pdery,pderz,xlap)
use field  ! (pot4,hpsi)
use grid   ! (nx,ny,nz,nxyz,dxyz)
use gridk  ! (px,py,pz)
use he4    ! (h2o2m4)
use rho    ! (*psi* & *psiv*)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)

implicit none

real    (kind=8), Allocatable :: hh(:)
real    (kind=8) :: deltat,fac
real    (kind=8) :: mu4,mumean
real    (kind=8) :: a1,sto,aux
complex (kind=8) :: auxc,ci

integer (kind=4) :: ix,iy,iz,is

If(.Not.Allocated(hh))Then
  Allocate(hh(3))
  hh(1)=hx  
  hh(2)=hy  
  hh(3)=hz  
EndIf


ci=(0.d0,1.d0)
fac=158.66624

Do is=1,3

  Call derivnD(1,nn,hh(is),is,psi,sto1c,Icon)
!
! now sto1c contains the is-component of Grad(Psi)
!

!    vel(is,:,:,:)=fac*dimag(sto1c/psi)
  Do iz=1, nz
    Do iy=1, ny
      Do ix=1, nx 
        aux =den(ix,iy,iz)
        If(aux.Gt.Min_rho.And.aux.Lt.Max_rho)Then
          vel(is,ix,iy,iz)=fac*dimag(sto1c(ix,iy,iz)/psi(ix,iy,iz))
        Else
          vel(is,ix,iy,iz)=0.d0
        EndIf  
      Enddo
    Enddo
  Enddo
Enddo


return

end subroutine velocity
!--------------------------------------------------------------------
!---                  Subroutine VForma                           ---
!--------------------------------------------------------------------

! Gives the FFT of the backflow term


subroutine vforma(nx,ny,nz,pmod,ulj)
implicit none

integer   (kind=4) :: nx,ny,nz
real      (kind=8) :: p
real      (kind=8) :: pmod(nx/2+1,ny,nz)
real      (kind=8) :: ulj(nx/2+1,ny,nz)

real      (kind=8)  :: pi,twopi
real      (kind=8)  :: aux1,aux2,aux3
integer   (kind=4)  :: ipx,ipy,ipz
real      (kind=8)  :: G11=-19.7544d0,G12=12.5616d0,G21=-0.2395d0,G22=0.0312d0,A1=1.023d0,A2=0.14912d0
real      (kind=8)  :: fac,pi3o2,pi2,qx
real      (kind=8)  :: B1,B2,f11,f12,f21,f22

pi    = 4.0d0*datan(1.0d0)
twopi = 2.0d0*pi
pi2=pi*pi
b1=pi2/a1
b2=pi2/a2
pi3o2=pi*Dsqrt(pi)
fac=pi3o2*0.5d0
f11=fac*(3.*g12+2.*g11*a1)/(a1**2*dsqrt(a1))
f12=-fac*2.*g12*pi2/(a1**3*dsqrt(a1))
!Open(unit=53,File='vforma.out')
!Write(53,'("f11,f12,b1..",1p,3E25.16)')f11,f12,b1
f21=fac*(3.*g22+2.*g21*a2)/(a2**2*dsqrt(a2))
f22=-fac*2.*g22*pi2/(a2**3*dsqrt(a2))
!Write(53,'("f21,f22,b2..",1p,3E25.16)')f21,f22,b2
!Call Flush(53)
!Close(53)

do ipz=1,nz
  do ipy=1,ny
    do ipx=1,nx/2+1
       p = pmod(ipx,ipy,ipz)
       qx=p
!      computed by ALberto, OK!
!       ulj(ipx,ipy,ipz) = (-7.187834473463007 - 637.5361220978009*qx**2)*dexp(-9.647707136939758*qx**2) &
!                   +(7.188703041045617 - 1339.0707262035585*qx**2)*dexp(-66.18565183133958*qx**2)
!      computed by Marti, OK!
       ulj(ipx,ipy,ipz) = (f11 + f12*qx**2)*dexp(-b1*qx**2) &
                   +(f21 + f22*qx**2)*dexp(-b2*qx**2)
    end do
  end do
end do
return

end subroutine vforma

