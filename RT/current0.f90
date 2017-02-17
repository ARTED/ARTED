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
!This file is "current.f90"
!This file contain one subroutine.
!SUBROUTINE current
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif

Subroutine current0_omp_KB
  use Global_Variables
  use timelog
#ifdef ARTED_USE_NVTX
  use nvtx
#endif
  implicit none
  real(8) :: jx,jy,jz

  NVTX_BEG('current_omp_KB()',2)
  call timelog_begin(LOG_CURRENT)
#ifndef _OPENACC
  call current0_omp_KB_impl(zu,jx,jy,jz)
#else
  call current0_acc_KB_impl(zu,jx,jy,jz)
#endif
  call current_result_impl(jx,jy,jz)
  call timelog_end(LOG_CURRENT)
  NVTX_END()
end Subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine current0_omp_KB_impl(zutmp,jxs,jys,jzs)
  use Global_Variables
  implicit none
  complex(8),intent(in) :: zutmp(0:NL-1,NBoccmax,NK_s:NK_e)
  real(8),intent(out)   :: jxs,jys,jzs

  integer :: ikb,ib,ik,i,j,ix,iy,iz,ia
  real(8) :: kr,jx,jy,jz,IaLxyz
  real(8) :: nabt(12)

  nabt( 1: 4) = nabx(1:4)
  nabt( 5: 8) = naby(1:4)
  nabt( 9:12) = nabz(1:4)

  IaLxyz = 1.0 / aLxyz

  jxs=0.d0
  jys=0.d0
  jzs=0.d0
!$omp parallel reduction(+:jxs,jys,jzs)

!Constructing nonlocal part
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
!$omp end do

!$omp do private(ikb,ik,ib,jx,jy,jz)
  do ikb=1,NKB
    ik=ik_table(ikb)
    ib=ib_table(ikb)

    call current_init_impl(ik,ib,zutmp(:,ib,ik),jx,jy,jz)
    call current_stencil_impl(ik,ib,zutmp(:,ib,ik),nabt,jx,jy,jz)
    call current_pseudo_impl(ik,ib,zutmp(:,ib,ik),IaLxyz,jx,jy,jz)

    jxs=jxs+jx
    jys=jys+jy
    jzs=jzs+jz
  enddo
!$omp end parallel
end Subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_LBLK
subroutine current0_acc_KB_impl(zutmp,jxs,jys,jzs)
  use Global_Variables
  use opt_variables
  implicit none
  complex(8),intent(in) :: zutmp(0:NL-1,NBoccmax,NK_s:NK_e)
  real(8),intent(out)   :: jxs,jys,jzs

  integer :: ikb,ib,ik,i,j,ix,iy,iz,ia
  real(8) :: kr,IaLxyz
  real(8) :: nabt(12)
  real(8) :: jx(NKB), jy(NKB), jz(NKB)
  integer :: ikb0, num_ikb1, ikb_s,ikb_e

  nabt( 1: 4) = nabx(1:4)
  nabt( 5: 8) = naby(1:4)
  nabt( 9:12) = nabz(1:4)

  IaLxyz = 1.0 / aLxyz

!$acc data pcopyin(zutmp) create(jx,jy,jz) pcopyout(ekr_omp) copyin(nabt) &
!$acc& pcopyin(jxyz,jxx,jyy,jzz,kAc,lx,ly,lz,Mps) 

!Constructing nonlocal part
!$acc kernels
!$acc loop collapse(2) independent gang
  do ik=NK_s,NK_e
  do ia=1,NI
!$acc loop independent vector(128)
  do j=1,Mps(ia)
    i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
    kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
    ekr_omp(j,ia,ik)=exp(zI*kr)
  end do
  end do
  end do
!$acc end kernels

  do ikb0 = 1, NKB, blk_nkb_current
    num_ikb1 = min(blk_nkb_current, NKB-ikb0+1)
    ikb_s = ikb0
    ikb_e = ikb0 + num_ikb1-1

    call current_init_impl_LBLK(zutmp(:,:,:), jx(:),jy(:),jz(:), ikb_s,ikb_e)
    call current_RT_stencil_impl_LBLK(jx(:),jy(:),jz(:), ikb_s,ikb_e)
    call current_pseudo_impl_LBLK(zutmp(:,:,:),IaLxyz,jx(:),jy(:),jz(:), ikb_s,ikb_e)
  enddo

!$acc kernels
  jxs=0.d0
  jys=0.d0
  jzs=0.d0
!$acc loop gang vector reduction(+:jxs,jys,jzs)
  do ikb = 1, NKB
    jxs=jxs+jx(ikb)
    jys=jys+jy(ikb)
    jzs=jzs+jz(ikb)
  enddo
!$acc end kernels

!$acc end data
end subroutine
#endif
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
