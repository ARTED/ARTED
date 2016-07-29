!
!  Copyright 2016 ARTED developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!This file is "dt_evolve.f90"
!This file contain one subroutine.
!Subroutine dt_evolve
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif

Subroutine dt_evolve_omp_KB(iter)
  use Global_Variables
  use timelog
#ifdef ARTED_USE_NVTX
  use nvtx
#endif
  use opt_variables
  implicit none
  integer    :: ik,ib,iter,ixyz
  integer    :: ia,j,i,ix,iy,iz
  real(8)    :: kr
  integer    :: thr_id,omp_get_thread_num,ikb

  NVTX_BEG('dt_evolve_omp_KB',1)
  call timelog_begin(LOG_DT_EVOLVE)

  NVTX_BEG('nonlocal part?',2)
  thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()

!Constructing nonlocal part ! sato
!$omp do private(ik,ia,j,i,ix,iy,iz,kr) collapse(2)
  do ik=NK_s,NK_e
  do ia=1,NI
  do j=1,Mps(ia)
    i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
    kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
    ekr_omp(j,ia,ik)=exp(zI*kr)
  end do
  end do
  end do
!$omp end parallel
  NVTX_END()

! yabana
  select case(functional)
  case('VS98','TPSS','TBmBJ')

!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu_GS(:,ib,ik)=zu(:,ib,ik)
  end do
    
  Vloc_t=Vloc

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass_t=tmass
    tjr_t=tjr
    tjr2_t=tjr2
  end if

  call dt_evolve_hpsi

  call psi_rho_RT
  call Hartree
  call Exc_Cor('RT')

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  Vloc=0.5d0*(Vloc+Vloc_t)

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass=0.5d0*(tmass+tmass_t)
    tjr=0.5d0*(tjr+tjr_t)
    tjr2=0.5d0*(tjr2+tjr2_t)
  end if

  if (Longi_Trans == 'Lo') then 
    call current_omp_KB
    javt(iter,:)=jav(:)
    Ac_ind(iter+1,:)=2*Ac_ind(iter,:)-Ac_ind(iter-1,:)-4*Pi*javt(iter,:)*dt**2
    if (Sym /= 1) then
      Ac_ind(iter+1,1)=0.d0
      Ac_ind(iter+1,2)=0.d0
    end if
    Ac_tot(iter+1,:)=Ac_ext(iter+1,:)+Ac_ind(iter+1,:)
  else if (Longi_Trans == 'Tr') then 
    Ac_tot(iter+1,:)=Ac_ext(iter+1,:)
  end if

!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu(:,ib,ik)=zu_GS(:,ib,ik)
  end do

  do ixyz=1,3
    kAc(:,ixyz)=kAc0(:,ixyz)+0.5*(Ac_tot(iter,ixyz)+Ac_tot(iter+1,ixyz))
  enddo



  end select
! yabana

!$acc data pcopy(zu)  
  NVTX_BEG('dt_evolve_hpsi',3)
  call dt_evolve_hpsi
  NVTX_END()

  NVTX_BEG('psi_rho_RT',4)
  call psi_rho_RT
  NVTX_END()
!$acc end data

  NVTX_BEG('Hartree',5)
  call Hartree
  NVTX_END()

! yabana
  NVTX_BEG('Exc_Cor',6)
  call Exc_Cor('RT')
  NVTX_END()
! yabana

  do ixyz=1,3
    kAc(:,ixyz)=kAc0(:,ixyz)+Ac_tot(iter,ixyz)
  enddo

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  call timelog_end(LOG_DT_EVOLVE)
  NVTX_END()

  return
End Subroutine dt_evolve_omp_KB
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine dt_evolve_omp_KB_MS
  use Global_Variables
  use timelog
  use opt_variables
  implicit none
  integer    :: ik,ib
  integer    :: ia,j,i,ix,iy,iz
  real(8)    :: kr
  integer    :: thr_id,omp_get_thread_num,ikb

  call timelog_begin(LOG_DT_EVOLVE)

  thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()

!Constructing nonlocal part ! sato
!$omp do private(ik,ia,j,i,ix,iy,iz,kr) collapse(2)
  do ik=NK_s,NK_e
  do ia=1,NI
  do j=1,Mps(ia)
    i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
    kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
    ekr_omp(j,ia,ik)=exp(zI*kr)
  end do
  end do
  end do
!$omp end parallel

! yabana
  select case(functional)
  case('VS98','TPSS','TBmBJ')


!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu_GS(:,ib,ik)=zu(:,ib,ik)
  end do

  Vloc_t=Vloc

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass_t=tmass
    tjr_t=tjr
    tjr2_t=tjr2
  end if

  call dt_evolve_hpsi

  call psi_rho_RT
  call Hartree
  call Exc_Cor('RT')

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  Vloc=0.5d0*(Vloc+Vloc_t)

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass=0.5d0*(tmass+tmass_t)
    tjr=0.5d0*(tjr+tjr_t)
    tjr2=0.5d0*(tjr2+tjr2_t)
  end if

!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu(:,ib,ik)=zu_GS(:,ib,ik)
  end do

  end select
! yabana

  call dt_evolve_hpsi

  call psi_rho_RT
  call Hartree
! yabana
  call Exc_Cor('RT')
! yabana

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  call timelog_end(LOG_DT_EVOLVE)

  return
End Subroutine dt_evolve_omp_KB_MS
