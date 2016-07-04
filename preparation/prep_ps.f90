!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of prep_ps.f90.
!
!  prep_ps.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  prep_ps.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with prep_ps.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!This file is "prep_ps.f90"
!This file conatain one soubroutine.
!SUBROUTINE prep_ps_periodic(property)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine prep_ps_periodic(property)
  use Global_Variables
  implicit none
  character(11) :: property
  integer :: ik,n,i,a,j,ix,iy,iz,lma,l,m,lm,ir,intr
  integer :: lma_tbl((Lmax+1)**2,NI)
  real(8) :: G2sq,s,Vpsl_l(NL),G2,Gd,Gr,x,y,z,r,ratio1,ratio2
  real(8) :: Ylm,dYlm,uVr(0:Lmax),duVr(0:Lmax)
  complex(8) :: Vion_G(NG_s:NG_e)
  !spline interpolation
  real(8) :: xx
  real(8) :: udVtbl_a(Nrmax,0:Lmax),dudVtbl_a(Nrmax,0:Lmax)
  real(8) :: udVtbl_b(Nrmax,0:Lmax),dudVtbl_b(Nrmax,0:Lmax)
  real(8) :: udVtbl_c(Nrmax,0:Lmax),dudVtbl_c(Nrmax,0:Lmax)
  real(8) :: udVtbl_d(Nrmax,0:Lmax),dudVtbl_d(Nrmax,0:Lmax)
  real(8),allocatable :: xn(:),yn(:),an(:),bn(:),cn(:),dn(:)  
  real(8) :: vloc_av
  if(property == 'not_initial') then
    deallocate(a_tbl,uV,duV,iuV,Jxyz,Jxx,Jyy,Jzz)
    deallocate(ekr) ! sato
    deallocate(ekr_omp)
  end if

! local potential
!$omp parallel
!$omp do private(ik,n,G2sq,s,r,i,vloc_av) collapse(2)
  do ik=1,NE
    do n=NG_s,NG_e
      G2sq=sqrt(Gx(n)**2+Gy(n)**2+Gz(n)**2)
      s=0.d0
      if (n == nGzero) then
        do i=2,NRloc(ik)
           r=0.5d0*(rad(i,ik)+rad(i-1,ik))
           vloc_av = 0.5d0*(vloctbl(i,ik)+vloctbl(i-1,ik))
           s=s+4*Pi*(r**2*vloc_av+r*Zps(ik))*(rad(i,ik)-rad(i-1,ik))
        enddo
      else
        do i=2,NRloc(ik)
          r=0.5d0*(rad(i,ik)+rad(i-1,ik))
          vloc_av = 0.5d0*(vloctbl(i,ik)+vloctbl(i-1,ik))
          s=s+4*Pi*sin(G2sq*r)/(G2sq)*(r*vloc_av+Zps(ik))*(rad(i,ik)-rad(i-1,ik))
        enddo
      endif
      dVloc_G(n,ik)=s
    enddo
  enddo
!$omp end do
!$omp end parallel

  rhoion_G=0.d0
!$omp parallel private(a)
  do a=1,NI
!$omp do private(n)
    do n=NG_s,NG_e
      rhoion_G(n)=rhoion_G(n)+Zps(Kion(a))/aLxyz*exp(-zI*(Gx(n)*Rion(1,a)+Gy(n)*Rion(2,a)+Gz(n)*Rion(3,a)))
    enddo
!$omp end do
  enddo
!$omp end parallel

  Vion_G=0.d0
!$omp parallel private(a,ik)
  do a=1,NI
    ik=Kion(a)
!$omp do private(n,G2,Gd)
    do n=NG_s,NG_e
      G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
      Gd=Gx(n)*Rion(1,a)+Gy(n)*Rion(2,a)+Gz(n)*Rion(3,a)
      Vion_G(n)=Vion_G(n)+dVloc_G(n,ik)*exp(-zI*Gd)/aLxyz
      if(n == nGzero) cycle
      Vion_G(n)=Vion_G(n)-4*Pi/G2*Zps(ik)*exp(-zI*Gd)/aLxyz
    enddo
!$omp end do
  enddo
!$omp end parallel

  Vpsl_l=0.d0
!$omp parallel private(n)
  do n=NG_s,NG_e
!$omp do private(i,Gr)
    do i=1,NL
      Gr=Gx(n)*Lx(i)*Hx+Gy(n)*Ly(i)*Hy+Gz(n)*Lz(i)*Hz
      Vpsl_l(i)=Vpsl_l(i)+Vion_G(n)*exp(zI*Gr)
    enddo
!$omp end do
  enddo
!$omp end parallel

  call MPI_ALLREDUCE(Vpsl_l,Vpsl,NL,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)

! nonlocal potential
  if (Myrank == 0 .and. property=='initial') then
    write(*,*) ''
    write(*,*) '============nonlocal grid data=============='
  endif
  do a=1,NI
    ik=Kion(a)
    j=0
    do ix=-2,2
    do iy=-2,2
    do iz=-2,2
      do i=1,NL
        x=Lx(i)*Hx-(Rion(1,a)+ix*aLx)
        y=Ly(i)*Hy-(Rion(2,a)+iy*aLy)
        z=Lz(i)*Hz-(Rion(3,a)+iz*aLz)
        r=sqrt(x*x+y*y+z*z)
        if (r<Rps(ik)) then
          j=j+1
        endif
      enddo
    enddo
    enddo
    enddo
    Mps(a)=j
    if (Myrank == 0 .and. property == 'initial') then
      write(*,*) 'a =',a,'Mps(a) =',Mps(a)
    endif
  end do
  Nps=maxval(Mps(:))
  allocate(Jxyz(Nps,NI),Jxx(Nps,NI),Jyy(Nps,NI),Jzz(Nps,NI))
  allocate(ekr(Nps,NI)) ! sato
  allocate(ekr_omp(Nps,NI,NK_s:NK_e))

  do a=1,NI
    ik=Kion(a)
    j=0
    do ix=-2,2
    do iy=-2,2
    do iz=-2,2
      do i=1,NL
        x=Lx(i)*Hx-(Rion(1,a)+ix*aLx)
        y=Ly(i)*Hy-(Rion(2,a)+iy*aLy)
        z=Lz(i)*Hz-(Rion(3,a)+iz*aLz)
        r=sqrt(x*x+y*y+z*z)
        if (r<Rps(ik)) then
          j=j+1
          if (j<=Nps) then
            Jxyz(j,a)=i
            Jxx(j,a)=ix
            Jyy(j,a)=iy
            Jzz(j,a)=iz
          endif
        endif
      enddo
    enddo
    enddo
    enddo
  end do

  lma=0
  do a=1,NI
    ik=Kion(a)
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lma=lma+1
      enddo
    enddo
  enddo
  Nlma=lma

  allocate(a_tbl(Nlma),uV(Nps,Nlma),iuV(Nlma),duV(Nps,Nlma,3))

  lma=0
  do a=1,NI
    ik=Kion(a)
    lm=0
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        lma=lma+1
        a_tbl(lma)=a
        lma_tbl(lm,a)=lma
      enddo
    enddo
  enddo

  do a=1,NI
    ik=Kion(a)
    do j=1,Mps(a)
      x=Lx(Jxyz(j,a))*Hx-(Rion(1,a)+Jxx(j,a)*aLx)
      y=Ly(Jxyz(j,a))*Hy-(Rion(2,a)+Jyy(j,a)*aLy)
      z=Lz(Jxyz(j,a))*Hz-(Rion(3,a)+Jzz(j,a)*aLz)
      r=sqrt(x**2+y**2+z**2)+1d-50
      do ir=1,NRps(ik)
        if(radnl(ir,ik).gt.r) exit
      enddo
      intr=ir-1
      if (intr.lt.0.or.intr.ge.NRps(ik))stop 'bad intr at prep_ps'
      ratio1=(r-radnl(intr,ik))/(radnl(intr+1,ik)-radnl(intr,ik))
      ratio2=1-ratio1
      do l=0,Mlps(ik)
        uVr(l)=ratio1*udVtbl(intr+1,l,ik)+ratio2*udVtbl(intr,l,ik)
        duVr(l)=ratio1*dudVtbl(intr+1,l,ik)+ratio2*dudVtbl(intr,l,ik)
      enddo
      lm=0
      do l=0,Mlps(ik)
        if(inorm(l,ik)==0) cycle
        do m=-l,l
          lm=lm+1
          uV(j,lma_tbl(lm,a))=uVr(l)*Ylm(x,y,z,l,m)
          duV(j,lma_tbl(lm,a),1)=duVr(l)*(x/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,1)
          duV(j,lma_tbl(lm,a),2)=duVr(l)*(y/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,2)
          duV(j,lma_tbl(lm,a),3)=duVr(l)*(z/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,3)
        enddo
      enddo
    enddo
    lm=0
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        iuV(lma_tbl(lm,a))=inorm(l,ik)
      enddo
    enddo
  enddo

  return
End Subroutine prep_ps_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
