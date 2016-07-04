!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of Hartree.f90.
!
!  Hartree.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Hartree.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with Hartree.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!This file is "Hartree.f90"
!This file contain one subroutine.
!SUBROUTINE Hartree
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Hartree
  use Global_Variables
  use timelog
  implicit none
  integer :: i,ix,iy,iz,n,nx,ny,nz
  real(8) :: G2

  call timelog_begin(LOG_HARTREE)

!$omp parallel 

!$omp do private(i)
  do i=1,NL
    rho_3D(Lx(i),Ly(i),Lz(i))=rho(i)
  end do
!$omp end do

!$omp do private(ix,iy,nz) collapse(3)
  do nz = -NLz/2,NLz/2-1
  do iy = 0,NLy-1
  do ix = 0,NLx-1
    f1(ix,iy,nz)=sum(eGzc(nz,:)*rho_3D(ix,iy,:))
  end do
  end do
  end do
!$omp end do

!$omp do private(ix,ny,nz) collapse(3)
  do nz = -NLz/2,NLz/2-1
  do ny = -NLy/2,NLy/2-1
  do ix = 0,NLx-1
    f2(ix,ny,nz)=sum(eGyc(ny,:)*f1(ix,:,nz))
  end do
  end do
  end do
!$omp end do

!$omp do private(nx,ny,nz) collapse(3)
  do nz = -NLz/2,NLz/2-1
  do ny = -NLy/2,NLy/2-1
  do nx = -NLx/2,NLx/2-1
    rhoe_G_3D(nx,ny,nz)=sum(eGxc(nx,:)*f2(:,ny,nz))/dble(NL)
  end do
  end do
  end do
!$omp end do

!$omp do private(nx,ny,nz,G2) collapse(3)
  do nz = -NLz/2,NLz/2-1
  do ny = -NLy/2,NLy/2-1
  do nx = -NLx/2,NLx/2-1
    rhoe_G_temp(nxyz(nx,ny,nz))=rhoe_G_3D(nx,ny,nz)
    G2=Gx(nxyz(nx,ny,nz))**2+Gy(nxyz(nx,ny,nz))**2+Gz(nxyz(nx,ny,nz))**2
    rhoe_G_3D(nx,ny,nz)=4*Pi/G2*rhoe_G_3D(nx,ny,nz)
  end do
  end do
  end do
!$omp end do

  rhoe_G_3D(0,0,0)=0.d0

!$omp do private(n)
  do n=NG_s,NG_e
    rhoe_G(n)=rhoe_G_temp(n)
  end do
!$omp end do

!$omp do private(nx,ny,iz) collapse(3)
  do iz = 0,NLz-1
  do ny = -NLy/2,NLy/2-1
  do nx = -NLx/2,NLx/2-1
    f3(nx,ny,iz)=sum(eGz(:,iz)*rhoe_G_3D(nx,ny,:))
  end do
  end do
  end do
!$omp end do

!$omp do private(nx,iy,iz) collapse(3)
  do iz = 0,NLz-1
  do iy = 0,NLy-1
  do nx = -NLx/2,NLx/2-1
    f4(nx,iy,iz)=sum(eGy(:,iy)*f3(nx,:,iz))
  end do
  end do
  end do
!$omp end do

!$omp do private(ix,iy,iz) collapse(3)
  do iz = 0,NLz-1
  do iy = 0,NLy-1
  do ix = 0,NLx-1
    Vh_3D(ix,iy,iz)=sum(eGx(:,ix)*f4(:,iy,iz))
  end do
  end do
  end do
!$omp end do

!$omp do private(i)
  do i=1,NL
    Vh(i)=Vh_3D(Lx(i),Ly(i),Lz(i))
  end do
!$omp end do

!$omp end parallel

  call timelog_end(LOG_HARTREE)

  return
End Subroutine Hartree
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Hartree_6
  use Global_Variables
  implicit none
  integer :: n,i
  real(8) :: Gr,G2,Vh_l(NL)

  Vh_l=0.d0
  do n=NG_s,NG_e
    G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
    rhoe_G(n)=0.d0
    do i=1,NL
      Gr=Gx(n)*Lx(i)*Hx+Gy(n)*Ly(i)*Hy+Gz(n)*Lz(i)*Hz
      rhoe_G(n)=rhoe_G(n)+rho(i)*exp(-zI*Gr)
    enddo
    rhoe_G(n)=rhoe_G(n)*Hxyz/aLxyz
    if(n == nGzero) cycle
    do i=1,NL
      Gr=Gx(n)*Lx(i)*Hx+Gy(n)*Ly(i)*Hy+Gz(n)*Lz(i)*Hz
      Vh_l(i)=Vh_l(i)+4*Pi/G2*real(rhoe_G(n)*exp(zI*Gr))
    enddo
  enddo
  call MPI_ALLREDUCE(Vh_l,Vh,NL,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)

  return
End Subroutine Hartree_6
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130

