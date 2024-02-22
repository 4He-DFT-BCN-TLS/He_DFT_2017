!-----------------------------------------------------------------------------
!--                          Subroutine evolo                              ---
!-----------------------------------------------------------------------------
!
subroutine evolo(deltat,mu4,mu4err,n4)

Use Para_DerivnD
use deriva ! (icon,npd,dxden,dyden,dzden,pderx,pdery,pderz,llap,xlap)
use vortex ! (to compute with a vortex)
use field  ! (pot4,hpsi,uext,limp)
use grid   ! (nx,ny,nz,nxyz,dxyz)
use gridk  ! (px,py,pz)
use he4    ! (h2o2m4)
use impur  !
use rho    ! (*psi*)
use util1  ! (vdt,nn,mmx,iw)
use work1  ! (temporal storage)

implicit none

integer (kind=4) :: n4
real    (kind=8) :: deltat,mu,muerr,n4real
real    (kind=8) :: mu4,mu4err,xcm,ycm,zcm,mu4err_gs
real    (kind=8), save :: mu4_gs

integer (kind=4) :: ix,iy,iz
real    (kind=8) :: hmean,a1
complex    (kind=8) :: caux,ci=(0.d0,1.0d0)
real    (kind=8) :: cnorm,sto,aux,dlaux,aux2,aux1


dlaux=dlog(2.0d0)


!.......................
!.. Laplacian of Psi ...
!.......................

!
!   icon   =  0 ! Take the derivative.
!   icon   = 20 ! Take the derivative. Use Wigner-Seitz conditions.
!   icon   =  8 ! Take the derivative. Use periodic conditions.
   
  Call derivnD(2,nn,hx,1,psi,sto1c,Icon)
  Call derivnD(2,nn,hy,2,psi,sto2c,Icon)
  Call derivnD(2,nn,hz,3,psi,sto3c,Icon)

  sto4c = -(sto1c+sto2c+sto3c)*h2o2m4 + pot4*psi
!
!     Here we will compute -Omega*|L_axis Psi>
!
  If(Lvortex.And.nv.Gt.1.And..Not.Lvortex_dipole)Then
    xcm=0.0d0
    ycm=0.0d0
    zcm=0.0d0
    Do ix=1, nx
      xcm = xcm + Sum(den(ix,:,:)*x(ix))
    EndDo 
    Do iy=1, ny
      ycm = ycm + Sum(den(:,iy,:)*y(iy))
    EndDo 
    Do iz=1, nz
      zcm = zcm + Sum(den(:,:,iz)*z(iz))
    EndDo 
    xcm = xcm*dxyz/n4
    ycm = ycm*dxyz/n4
    zcm = zcm*dxyz/n4
    caux = (0.d0, 0.d0)      
    If(Vortex_axis.Eq.'Z')Then
      Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
      Call derivnD(1,nn,hy,2,psi,sto2c,Icon)
      Do iz=1, nz
        Do iy=1, ny
          Do ix=1, nx
            sto5c(ix,iy,iz) = Ci*((y(iy)-ycm)*sto1c(ix,iy,iz) - (x(ix)-xcm)*sto2c(ix,iy,iz)) 
          EndDo
        EndDo
      EndDo
    Endif
    If(Vortex_axis.Eq.'Y')Then
      Call derivnD(1,nn,hx,1,psi,sto1c,Icon)
      Call derivnD(1,nn,hz,3,psi,sto3c,Icon)
      Do iz=1, nz
        Do iy=1, ny
          Do ix=1, nx
            sto5c(ix,iy,iz) = Ci*((x(ix)-xcm)*sto3c(ix,iy,iz) - (z(iz)-zcm)*sto1c(ix,iy,iz)) 
          EndDo
        EndDo
      EndDo
    Endif
    If(Vortex_axis.Eq.'X')Then
      Call derivnD(1,nn,hy,2,psi,sto2c,Icon)
      Call derivnD(1,nn,hz,3,psi,sto3c,Icon)
      Do iz=1, nz
        Do iy=1, ny
          Do ix=1, nx
            sto5c(ix,iy,iz) = Ci*((z(iz)-zcm)*sto2c(ix,iy,iz) - (y(iy)-ycm)*sto3c(ix,iy,iz)) 
          EndDo
        EndDo
      EndDo
    Endif
    sto4c = sto4c -omega*sto5c
  Endif        
!.................................................................. hpsi  = (T+U)*Psi_4 (With impurity)


   hpsi = sto4c
   hmean  = sum(conjg(psi)*hpsi)/sum(den)    ! Average over all the mu4's
   mu4err = abs(1.0d0-mu4/hmean)             ! Relative error
  
!   Write(6,'("From evolo: hmean, mu4, mu4err....",1p,3E15.6)')hmean,mu4,mu4err

   a1=1.0d0-deltat*(hmean-mu4)
   a1=deltat/a1

   If(.Not.Lbulk)Then
     mu4    = hmean
   Endif

   hpsi =a1*(hpsi-mu4*psi)

   if(iron.ne.0) then
     do ix=1,iron
       call ironingc(hpsi,nx,ny,nz)       ! Smoothing  (H-mu)*Psi
     end do
   end if

     psinw = psi - hpsi               &
           + vdt(1) * (psi  - psi1)   &
           + vdt(2) * (psi1 - psi2)
     psi2  = psi1
     psi1  = psi

   hpsi = sto4c

   If(Lbulk)Then
     cnorm=1.0d0
     n4real=sum(den)*dxyz
     n4=n4real+0.5
   Else
     cnorm = sqrt(n4/(sum(conjg(psinw)*psinw)*dxyz))
   Endif

     psi = psinw*cnorm
     den = conjg(psi)*psi

   If(LOrtogonal)Then
     sto5c = hpsi

     Call derivnD(2,nn,hx,1,psi_gs,sto1c,Icon)
     Call derivnD(2,nn,hy,2,psi_gs,sto2c,Icon)
     Call derivnD(2,nn,hz,3,psi_gs,sto3c,Icon)

     sto4c = -(sto1c+sto2c+sto3c)*h2o2m4 + pot4*psi_gs

     hpsi = sto4c
     hmean  = sum(conjg(psi_gs)*hpsi)/sum(conjg(psi_gs)*psi_gs)    ! Average over all the mu4's
     mu4err_gs = abs(1.0d0-mu4_gs/hmean)             ! Relative error
  
!     Write(6,'("From evolo: hmean, mu4_gs, mu4err_gs....",1p,3E15.6)')hmean,mu4_gs,mu4err_gs

     a1=1.0d0-deltat*(hmean-mu4_gs)
     a1=deltat/a1
     mu4_gs    = hmean

     hpsi =a1*(hpsi-mu4_gs*psi_gs)

     if(iron.ne.0) then
       do ix=1,iron
         call ironingc(hpsi,nx,ny,nz)       ! Smoothing  (H-mu)*Psi
       end do
     end if

       psi_gsnw = psi_gs -  hpsi                 &
                + vdt(1) * (psi_gs  - psi_gs1)   &
                + vdt(2) * (psi_gs1 - psi_gs2)
       psi_gs2  = psi_gs1
       psi_gs1  = psi_gs
!       
!  Aqui normalizamos la Psi del g.s. a 1   
!
     cnorm  = Sqrt(1.d0/(Sum(Conjg(psi_gsnw)*psi_gsnw)*dxyz))
     psi_gs = psi_gsnw*cnorm
!     
!  Aqui ortogonalizamos la Psi con la del g.s.     
!     
     caux  = Sum(Conjg(Psi_gs)*Psi)*dxyz
!     Write(6,'("From evolo: Crossing integral, <Psi|Psi_gs>....:",1p,2E18.10)')caux
     Psi   = Psi - Caux*Psi_gs
     Cnorm = Sqrt(n4/(Sum(Conjg(Psi)*Psi)*dxyz))
     psi    = psi*Cnorm
     den    = Conjg(psi)*psi

     hpsi = sto5c
   Endif
   If(L_cilindrical_average)Then
     Call den_cil_average(n4)
   EndIf        
return
end
