!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of Gram_Schmidt.f90.
!
!  Gram_Schmidt.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Gram_Schmidt.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with Gram_Schmidt.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!This file is "Gram_Schmidt.f90"
!This file contain one subroutine.
!Subroutine Gram_Schmidt
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Gram_Schmidt
  use Global_Variables
  use timelog
  implicit none
  integer :: ik,ib,ibt
  real(8) :: s
  complex(8) :: zov

  call timelog_begin(LOG_GRAM_SCHMIDT)
!$omp parallel do private(ib,ibt,zov,s)
  do ik=NK_s,NK_e
  do ib=1,NB
    do ibt=1,ib-1
      zov=sum(conjg(zu_GS(:,ibt,ik))*zu_GS(:,ib,ik))*Hxyz
      zu_GS(:,ib,ik)=zu_GS(:,ib,ik)-zu_GS(:,ibt,ik)*zov
    enddo
    s=sum(abs(zu_GS(:,ib,ik))**2)*Hxyz
    zu_GS(:,ib,ik)=zu_GS(:,ib,ik)/sqrt(s)
  enddo
  enddo
  call timelog_end(LOG_GRAM_SCHMIDT)

  return
End Subroutine Gram_Schmidt
