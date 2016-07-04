!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of Fourier_tr.f90.
!
!  Fourier_tr.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Fourier_tr.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with Fourier_tr.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Fourier_tr
  use Global_Variables
  implicit none
  integer :: ihw,ixyz,iter
  real(8) :: hw,tt
  complex(8) :: jav_w(3),E_ext_w(3),E_tot_w(3)
  complex(8) :: zsigma_w(3),zeps(3)

  if (Myrank == 0) then
    open(7,file=file_epse)
  endif

  do ihw=1,Nomega
    hw=ihw*domega
    jav_w=0.d0; E_ext_w=0.d0; E_tot_w=0.d0
    do iter=0,Nt
      tt=(iter+0.5)*dt
      jav_w(:)=jav_w(:)+javt(iter,:)*exp(zI*hw*tt)*smoothing_t(tt)
      E_ext_w(:)=E_ext_w(:)+E_ext(iter,:)*exp(zI*hw*tt)*smoothing_t(tt)
      E_tot_w(:)=E_tot_w(:)+E_tot(iter,:)*exp(zI*hw*tt)*smoothing_t(tt)
    enddo
    jav_w(:)=jav_w(:)*dt; E_ext_w(:)=E_ext_w(:)*dt; E_tot_w(:)=E_tot_w(:)*dt
    if (ext_field == 'LR') then
      if (Longi_Trans == 'Lo')  then 
        zeps(:)=1.d0/(1.d0-E_tot_w(:)/dAc)
      else if (Longi_Trans == 'Tr') then
        zsigma_w(:)=jav_w(:)/dAc
        zeps=1.d0+zI*4.d0*pi*zsigma_w(:)/hw
      end if
    end if
      
    if (Myrank == 0) then
      if (ext_field == 'LR' .and. Longi_Trans == 'Lo') then
        write(7,'(1x,f13.7,6f22.14)') hw&
             &,(real(zeps(ixyz)),ixyz=1,3)&
             &,(imag(zeps(ixyz)),ixyz=1,3)
      else if (ext_field == 'LR' .and. Longi_Trans == 'Tr') then
        write(7,'(1x,f13.7,12f22.14)') hw&
             &,(real(zsigma_w(ixyz)),ixyz=1,3)&
             &,(imag(zsigma_w(ixyz)),ixyz=1,3)&
             &,(real(zeps(ixyz)),ixyz=1,3)&
             &,(imag(zeps(ixyz)),ixyz=1,3)
      else
        write(7,'(1x,f13.7,18f22.14)') hw&
             &,(real(jav_w(ixyz)),ixyz=1,3)&
             &,(imag(jav_w(ixyz)),ixyz=1,3)&
             &,(real(E_ext_w(ixyz)),ixyz=1,3)&
             &,(imag(E_ext_w(ixyz)),ixyz=1,3)&
             &,(real(E_tot_w(ixyz)),ixyz=1,3)&
             &,(imag(E_tot_w(ixyz)),ixyz=1,3)
      endif
    endif
  enddo
 
  if (Myrank == 0) then
    close(7)
  endif

  return
Contains
!======
!======
Function smoothing_t(tt)
  implicit none
  real(8),intent(IN) :: tt
  real(8) :: smoothing_t
  
  smoothing_t=1.d0-3.d0*(tt/(Nt*dt))**2+2.d0*(tt/(Nt*dt))**3

  return
End Function smoothing_t
!======
!======
End Subroutine Fourier_tr
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
