!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of fd_coef.f90.
!
!  fd_coef.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  fd_coef.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with fd_coef.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!This file is "fd_coef.f90"
!This file contain one subroutine.
!SUBROUTINE fd_coef
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine fd_coef
  use Global_Variables
  implicit none
  integer :: i,j,k,n
  real(8) :: t,s,s0

  lap=0.d0
  nab=0.d0

  n=Nd
  s=0.d0
  do i=-n,n
    if ( i==0 ) cycle
    do j=-n,n
      if ( j==i .or. j==0 ) cycle
      s=s+1.d0/dble(i*j)
    end do
  end do
  lap(0)=s

  do j=1,n
    t=1.d0
    do i=-n,n
      if ( i==j ) cycle
      t=t*(j-i)
    end do
    s=0.d0
    do k=-n,n
      if ( k==j .or. k==0 ) cycle
      s0=1.d0
      do i=-n,n
        if ( i==k .or. i==j .or. i==0 ) cycle
        s0=s0*(-i)
      end do
      s=s+s0
    end do
    lap( j)=2.d0*s/t
    lap(-j)=lap(j)
  end do

  do j=1,n
    t=1.d0
    do i=-n,n
      if ( i==j ) cycle
      t=t*(j-i)
    end do
    s=1.d0
    do i=-n,n
      if ( i==j .or. i==0 ) cycle
      s=s*(-i)
    end do
    nab( j)=s/t
    nab(-j)=-nab(j)
  end do

  lapx(-Nd:Nd)=lap(-Nd:Nd)/Hx**2
  lapy(-Nd:Nd)=lap(-Nd:Nd)/Hy**2
  lapz(-Nd:Nd)=lap(-Nd:Nd)/Hz**2
  nabx(-Nd:Nd)=nab(-Nd:Nd)/Hx
  naby(-Nd:Nd)=nab(-Nd:Nd)/Hy
  nabz(-Nd:Nd)=nab(-Nd:Nd)/Hz

  return
End Subroutine fd_coef
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
