Subroutine den_cil_average(n4)
Use rho
Use grid
Use field
Use vortex
Implicit none
Integer  (kind=4) :: n4,ix,iy,iz,ir
Real     (kind=8) :: x2,y2,z2,r,aux
Real     (kind=8) :: xcm,ycm,zcm
denr  = 0.d0
normr = 0
!
!  We compute the C.M. position
!
xcm = 0.d0
ycm = 0.d0
zcm = 0.d0
Do ix=1,nx
  xcm = xcm + sum(den(ix,:,:)*x(ix))
EndDo
Do iy=1,ny
  ycm = ycm + sum(den(:,iy,:)*y(iy))
EndDo  
Do iz=1,nz
  zcm = zcm + sum(den(:,:,iz)*z(iz))
EndDo  
xcm=xcm*dxyz/n4
ycm=ycm*dxyz/n4
zcm=zcm*dxyz/n4
!
! We perform a cilindrial average
!
If(Vortex_axis.Eq.'Z')Then      
  Do iz=1,nz
    Do iy=1,ny
      y2=(y(iy)-ycm)**2
      Do ix=1,nx
        x2=(x(ix)-xcm)**2
        r=Dsqrt(x2+y2)
        ir= (r/hr + 1.5)
        If(ir.Gt.nr)Cycle
        denr(ir,iz)=denr(ir,iz) + den(ix,iy,iz)
        normr(ir,iz)=normr(ir,iz)+1
      Enddo
    Enddo
    Do ir=1,nr
      if(normr(ir,iz).Ne.0)denr(ir,iz)=denr(ir,iz)/normr(ir,iz)
    Enddo
    Call Spls3(rr,denr(:,iz),nr,r,aux,1,Adenr(:,iz),Bdenr(:,iz),0,0)
  Enddo
  Do iz=1,nz
    Do iy=1,ny
      y2=(y(iy)-ycm)**2
      Do ix=1,nx
        x2=(x(ix)-xcm)**2
        r=Dsqrt(x2+y2)
        If(r.Le.rr(nr))Then
          Call Spls3(rr,denr(:,iz),nr,r,aux,1,Adenr(:,iz),Bdenr(:,iz),1,1)
          den(ix,iy,iz) = aux
        Else
          den(ix,iy,iz) = 0.d0
        EndIf        
      Enddo
    Enddo
  Enddo
ElseIf(Vortex_axis.Eq.'Y')Then  
  Do iy=1,ny
    Do iz=1,nz
      z2=(z(iz)-zcm)**2
      Do ix=1,nx
        x2=(x(ix)-xcm)**2
        r=Dsqrt(x2+z2)
        ir= (r/hr + 1.5)
        If(ir.Gt.nr)Cycle
        denr(ir,iy)=denr(ir,iy) + den(ix,iy,iz)
        normr(ir,iy)=normr(ir,iy)+1
      Enddo
    Enddo
    Do ir=1,nr
      if(normr(ir,iy).Ne.0)denr(ir,iy)=denr(ir,iy)/normr(ir,iy)
    Enddo
    Call Spls3(rr,denr(:,iy),nr,r,aux,1,Adenr(:,iy),Bdenr(:,iy),0,0)
  Enddo
  Do iy=1,ny
    Do iz=1,ny
      z2=(z(iz)-zcm)**2
      Do ix=1,nx
        x2=(x(ix)-xcm)**2
        r=Dsqrt(x2+z2)
        If(r.Le.rr(nr))Then
          Call Spls3(rr,denr(:,iy),nr,r,aux,1,Adenr(:,iy),Bdenr(:,iy),1,1)
          den(ix,iy,iz) = aux
        Else
          den(ix,iy,iz) = 0.d0
        EndIf        
      Enddo
    Enddo
  Enddo
Else
  Do ix=1,nx
    Do iz=1,nz
      z2=(z(iz)-zcm)**2
      Do iy=1,ny
        y2=(y(iy)-ycm)**2
        r=Dsqrt(y2+z2)
        ir= (r/hr + 1.5)
        If(ir.Gt.nr)Cycle
        denr(ir,ix)=denr(ir,ix) + den(ix,iy,iz)
        normr(ir,ix)=normr(ir,ix)+1
      Enddo
    Enddo
    Do ir=1,nr
      if(normr(ir,ix).Ne.0)denr(ir,ix)=denr(ir,ix)/normr(ir,ix)
    Enddo
    Call Spls3(rr,denr(:,ix),nr,r,aux,1,Adenr(:,ix),Bdenr(:,ix),0,0)
  Enddo
  Do ix=1,nx
    Do iz=1,ny
      z2=(z(iz)-zcm)**2
      Do iy=1,ny
        y2=(y(iy)-ycm)**2
        r=Dsqrt(y2+z2)
        If(r.Le.rr(nr))Then
          Call Spls3(rr,denr(:,ix),nr,r,aux,1,Adenr(:,ix),Bdenr(:,ix),1,1)
          den(ix,iy,iz) = aux
        Else
          den(ix,iy,iz) = 0.d0
        EndIf        
      Enddo
    Enddo
  Enddo
EndiF  
aux = n4/(sum(den)*dxyz)
den = Max((den*aux),1.d-99)
Return
End
