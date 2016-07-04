!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of Fermi_Dirac_distribution.f90.
!
!  Fermi_Dirac_distribution.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Fermi_Dirac_distribution.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with Fermi_Dirac_distribution.f90.  If not, see <http://www.gnu.org/licenses/>.
!
Subroutine Fermi_Dirac_distribution
  use Global_Variables
  implicit none
  real(8) :: chemical_potential
  real(8) :: chem_max,chem_min,elec_num
  real(8) :: beta_FD
  real(8),allocatable :: occ_l(:,:),esp_l(:,:),esp_temp(:,:)
  real(8) s,st
  integer :: ik,ib
  real(8) :: timer1,timer2
  
  timer1=MPI_WTIME()
  beta_FD=1d0/(KbTev/(2d0*Ry))
  allocate(occ_l(NB,NK),esp_l(NB,NK),esp_temp(NB,NK))
  occ_l=0d0 ; esp_l=0d0
  esp_l(:,NK_s:NK_e)=esp(:,NK_s:NK_e)
  call MPI_ALLREDUCE(esp_l,esp_temp,NB*NK,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)
  chem_max=maxval(esp_temp)
  chem_min=minval(esp_temp)
  if(myrank == 0)then
    write(*,*)'max esp =',chem_max
    write(*,*)'min esp =',chem_min
  end if
  chemical_potential=0.5d0*(chem_max+chem_min)
  do 

    do ik=NK_s,NK_e
      do ib=1,NB
        occ_l(ib,ik)=2.d0/(NKx*NKy*NKz)*wk(ik) &
          &/(exp(beta_FD*(esp(ib,ik)-chemical_potential))+1d0)
      end do
    end do

    st=sum(occ_l(:,NK_s:NK_e))
    call MPI_ALLREDUCE(st,s,1,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)
    elec_num=s

    if(abs(elec_num-dble(Nelec)) <= 1d-6)exit
    if(elec_num-dble(Nelec) > 0d0)then
      chem_max=chemical_potential
      chemical_potential=0.5d0*(chem_max+chem_min)
    else
      chem_min=chemical_potential
      chemical_potential=0.5d0*(chem_max+chem_min)
    end if
    if(chem_max-chem_min < 1d-10)exit

  end do

  call MPI_ALLREDUCE(occ_l,occ,NB*NK,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)
  st=sum(occ_l(Nelec/2+1:NB,NK_s:NK_e))
  call MPI_ALLREDUCE(st,s,1,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)

  timer2=MPI_WTIME()
  if(myrank == 0)then
    write(*,*)'Fermi-Dirac dist. time=',timer2-timer1,'sec'
    write(*,*)'chemical potential =',chemical_potential
    write(*,*)'elec_num =',elec_num
    write(*,*)'excited electron =',s
!    open(126,file='occ.dat')
!    do ik=1,NK
!      do ib=1,NB
!        write(126,*)ik,ib,occ(ib,ik)
!      end do
!    end do
!    close(126)
  end if

  return
end Subroutine Fermi_Dirac_distribution
