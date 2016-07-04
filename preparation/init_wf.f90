!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of init_wf.f90.
!
!  init_wf.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  init_wf.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with init_wf.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!This file is "init_wf.f90"
!This file contain two subroutines
!SUBROUTINE init_wf
!SUBROUTINE quickrnd(iseed,rnd)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine init_wf
  use Global_Variables
  implicit none
  integer :: iseed,ib,ik,i
  real(8) :: r2,x1,y1,z1,rnd

  zu_GS=0.d0
  iseed=123
  do ik=NK_s,NK_e
    do ib=1,NB
      call quickrnd(iseed,rnd)
      x1=aLx*rnd
      call quickrnd(iseed,rnd)
      y1=aLy*rnd
      call quickrnd(iseed,rnd)
      z1=aLz*rnd
      do i=1,NL
        r2=(Lx(i)*Hx-x1)**2+(Ly(i)*Hy-y1)**2+(Lz(i)*Hz-z1)**2
        zu_GS(i,ib,ik)=exp(-0.5d0*r2)
      enddo
    enddo
  enddo

  return
End Subroutine init_wf
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine quickrnd(iseed,rnd)
  implicit none
  integer,parameter :: im=6075,ia=106,ic=1283
  integer :: iseed
  real(8) :: rnd

  iseed=mod(iseed*ia+ic,im)
  rnd=dfloat(iseed)/dfloat(im)

  return
End Subroutine quickrnd
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130

