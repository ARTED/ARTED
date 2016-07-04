!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of sp_energy.f90.
!
!  sp_energy.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  sp_energy.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with sp_energy.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!This file is "sp_energy.f90"
!This file contain one subroutine.
!SUBROUTINE sp_energy
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine sp_energy
  use Global_Variables
  implicit none
  integer :: ik,ib
  real(8) :: esp_l(NB,NK)
! sato
  integer :: ia,j,i,ix,iy,iz
  real(8) :: kr

  esp_l=0.d0
  do ik=NK_s,NK_e
!Constructing nonlocal part ! sato
  do ia=1,NI
    do j=1,Mps(ia)
      i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
      kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
      ekr(j,ia)=exp(zI*kr)
    enddo
  enddo
    do ib=1,NB
      tpsi(1:NL)=zu_GS(1:NL,ib,ik)
      call hpsi(ik)
      esp_l(ib,ik)=sum(conjg(zu_GS(:,ib,ik))*htpsi)*Hxyz
    enddo
  enddo

  CALL MPI_ALLREDUCE(esp_l,esp,NB*NK,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)

  return
End Subroutine sp_energy
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine sp_energy_omp
  use Global_Variables
  use timelog
  implicit none
  integer :: ik,ib
  real(8) :: esp_l(NB,NK)
! sato
  integer :: ia,j,i,ix,iy,iz
  real(8) :: kr
! omp
  integer :: thr_id,omp_get_thread_num

  call timelog_begin(LOG_SP_ENERGY)
  esp_l=0.d0
  thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ik,ia,j,i,ix,iy,iz,kr,ib)
  do ik=NK_s,NK_e
!Constructing nonlocal part ! sato
  do ia=1,NI
    do j=1,Mps(ia)
      i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
      kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
      ekr_omp(j,ia,ik)=exp(zI*kr)
    enddo
  enddo
    do ib=1,NB
      tpsi_omp(1:NL,thr_id)=zu_GS(1:NL,ib,ik)
      call hpsi_omp_KB(ik,tpsi_omp(:,thr_id),ttpsi_omp(:,thr_id),htpsi_omp(:,thr_id))
      esp_l(ib,ik)=sum(conjg(zu_GS(:,ib,ik))*htpsi_omp(:,thr_id))*Hxyz
    enddo
  enddo

!$omp end parallel

  call timelog_begin(LOG_ALLREDUCE)
  CALL MPI_ALLREDUCE(esp_l,esp,NB*NK,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)
  call timelog_end(LOG_ALLREDUCE)

  call timelog_end(LOG_SP_ENERGY)

  return
End Subroutine sp_energy_omp
