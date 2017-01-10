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


#ifdef ARTED_CURRENT_OPTIMIZED
subroutine current_omp_KB_ST(ib,ik,A)
  use Global_Variables
  use opt_variables
  implicit none
  integer,intent(in)    :: ib,ik
  complex(8),intent(in) :: A(0:NL-1)
  real(8) :: nabt(12),zx,zy,zz

  nabt( 1: 4) = nabx(1:4)
  nabt( 5: 8) = naby(1:4)
  nabt( 9:12) = nabz(1:4)

  zx = 0.d0
  zy = 0.d0
  zz = 0.d0
  call current_stencil(nabt,A,zx,zy,zz)

  zcx(ib,ik)=zx
  zcy(ib,ik)=zy
  zcz(ib,ik)=zz
end subroutine
#endif
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_CURRENT_OPTIMIZED
#ifdef ARTED_LBLK
subroutine current_acc_KB_ST_LBLK(A, ikb_s,ikb_e)
  use Global_Variables
  use opt_variables
  implicit none
  complex(8),intent(in) :: A(0:NL-1, NBoccmax, NK_s:NK_e)
  integer,intent(in)    :: ikb_s,ikb_e

  call current_stencil_LBLK(A(:,:,:), ikb_s,ikb_e)
end subroutine
#endif
#endif
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine current_GS_omp_KB
  use Global_Variables
  use timelog
  implicit none
  real(8) :: jx,jy,jz

  call timelog_begin(LOG_CURRENT)
  call current_GS_omp_KB_impl(zu_GS,jx,jy,jz)
  call timelog_end(LOG_CURRENT)
  call current_result_impl(jx,jy,jz)
end Subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine current_omp_KB
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
  call current_omp_KB_impl(zu,jx,jy,jz)
#else
  call current_acc_KB_impl(zu,jx,jy,jz)
#endif
  call timelog_end(LOG_CURRENT)
  call current_result_impl(jx,jy,jz)
  NVTX_END()
end Subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine current_init_impl(ik,ib,zutmp,jx,jy,jz)
  use Global_Variables
  implicit none
  integer,intent(in)    :: ik,ib
  complex(8),intent(in) :: zutmp(0:NL-1)
  real(8),intent(out)   :: jx,jy,jz

  integer :: i
  real(8) :: jt

  jt=0.d0
!dir$ vector aligned
  do i=0,NL-1
    jt=jt+real(zutmp(i))**2+imag(zutmp(i))**2
  end do

  jx=occ(ib,ik)*kAc(ik,1)*jt
  jy=occ(ib,ik)*kAc(ik,2)*jt
  jz=occ(ib,ik)*kAc(ik,3)*jt
end Subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_LBLK
Subroutine current_init_impl_LBLK(zutmp,jx,jy,jz, ikb_s,ikb_e)
  use Global_Variables
  implicit none
  complex(8),intent(in) :: zutmp(0:NL-1, NBoccmax,NK_s:NK_e)
  real(8),intent(out)   :: jx(NKB),jy(NKB),jz(NKB)
  integer,intent(in)    :: ikb_s,ikb_e

  integer :: i
  real(8) :: jt
  integer :: ikb, ik,ib
!$acc data pcopy(jx,jy,jz) pcopyin(zutmp,occ,kAc, ik_table,ib_table)
!$acc kernels
!$acc loop gang
  do ikb = ikb_s, ikb_e
    ik=ik_table(ikb)
    ib=ib_table(ikb)

    jt=0.d0
!$acc loop vector(256) reduction(+:jt)
    do i=0,NL-1
      jt=jt+real(zutmp(i,ib,ik))**2+imag(zutmp(i,ib,ik))**2
    end do

    jx(ikb)=occ(ib,ik)*kAc(ik,1)*jt
    jy(ikb)=occ(ib,ik)*kAc(ik,2)*jt
    jz(ikb)=occ(ib,ik)*kAc(ik,3)*jt
  enddo
!$acc end kernels
!$acc end data
end Subroutine
#endif ! ARTED_LBLK
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine current_stencil_impl(ik,ib,zutmp,nabt,jx,jy,jz)
  use Global_Variables
  use opt_variables
  implicit none
  integer,intent(in)    :: ik,ib
  complex(8),intent(in) :: zutmp(0:NL-1)
  real(8),intent(in)    :: nabt(12)
  real(8),intent(inout) :: jx,jy,jz
  real(8) :: jxt,jyt,jzt

  jxt=0.d0
  jyt=0.d0
  jzt=0.d0
  call current_stencil(nabt,zutmp,jxt,jyt,jzt)

  jx=jx+jxt*occ(ib,ik)
  jy=jy+jyt*occ(ib,ik)
  jz=jz+jzt*occ(ib,ik)
end Subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_CURRENT_OPTIMIZED
Subroutine current_RT_stencil_impl(ik,ib,jx,jy,jz)
  use Global_Variables
  use opt_variables
  implicit none
  integer,intent(in)    :: ik,ib
  real(8),intent(inout) :: jx,jy,jz

  jx=jx+zcx(ib,ik)*occ(ib,ik)
  jy=jy+zcy(ib,ik)*occ(ib,ik)
  jz=jz+zcz(ib,ik)*occ(ib,ik)
end Subroutine
#endif
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_LBLK
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine current_RT_stencil_impl_LBLK(jx,jy,jz, ikb_s,ikb_e)
  use Global_Variables
  use opt_variables
  implicit none
  real(8),intent(inout) :: jx(NKB),jy(NKB),jz(NKB)
  integer,intent(in)    :: ikb_s,ikb_e

  integer    :: ikb, ik,ib
!$acc data pcopy(jx,jy,jz) pcopyin(zcx,zcy,zcz,occ, ik_table,ib_table)
!$acc kernels
  do ikb = ikb_s, ikb_e
    ik=ik_table(ikb)
    ib=ib_table(ikb)

    jx(ikb)=jx(ikb)+zcx(ib,ik)*occ(ib,ik)
    jy(ikb)=jy(ikb)+zcy(ib,ik)*occ(ib,ik)
    jz(ikb)=jz(ikb)+zcz(ib,ik)*occ(ib,ik)
  enddo
!$acc end kernels
!$acc end data
end Subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#endif
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine current_pseudo_impl(ik,ib,zutmp,IaLxyz,jx,jy,jz)
  use Global_Variables
  use omp_lib
  use opt_variables
  use timelog
  implicit none
  integer,intent(in)    :: ik,ib
  complex(8),intent(in) :: zutmp(NL)
  real(8),intent(in)    :: IaLxyz
  real(8),intent(out)   :: jx,jy,jz

  integer    :: ilma,ia,j,i,ix,iy,iz
  real(8)    :: x,y,z
  real(8)    :: jxt,jyt,jzt
  complex(8) :: uVpsi,uVpsix,uVpsiy,uVpsiz

  jxt=0d0
  jyt=0d0
  jzt=0d0

  do ilma=1,Nlma
    ia=a_tbl(ilma)
    uVpsi=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
    do j=1,Mps(ia)
      i=Jxyz(j,ia)

      ix=Jxx(j,ia); x=Lx(i)*Hx-ix*aLx
      iy=Jyy(j,ia); y=Ly(i)*Hy-iy*aLy
      iz=Jzz(j,ia); z=Lz(i)*Hz-iz*aLz

      uVpsi =uVpsi +uV(j,ilma)*ekr_omp(j,ia,ik)  *zutmp(i)
      uVpsix=uVpsix+uV(j,ilma)*ekr_omp(j,ia,ik)*x*zutmp(i)
      uVpsiy=uVpsiy+uV(j,ilma)*ekr_omp(j,ia,ik)*y*zutmp(i)
      uVpsiz=uVpsiz+uV(j,ilma)*ekr_omp(j,ia,ik)*z*zutmp(i)
    end do
    uVpsi =uVpsi *Hxyz*iuV(ilma)
    uVpsix=uVpsix*Hxyz
    uVpsiy=uVpsiy*Hxyz
    uVpsiz=uVpsiz*Hxyz
    jxt=jxt+occ(ib,ik)*IaLxyz*2*imag(conjg(uVpsix)*uVpsi)
    jyt=jyt+occ(ib,ik)*IaLxyz*2*imag(conjg(uVpsiy)*uVpsi)
    jzt=jzt+occ(ib,ik)*IaLxyz*2*imag(conjg(uVpsiz)*uVpsi)
  enddo

  jx=jx*Hxyz*IaLxyz+jxt
  jy=jy*Hxyz*IaLxyz+jyt
  jz=jz*Hxyz*IaLxyz+jzt
end Subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_LBLK
Subroutine current_pseudo_impl_LBLK(zutmp,IaLxyz,jx,jy,jz, ikb_s,ikb_e)
  use Global_Variables
  use omp_lib
  use opt_variables
  use timelog
  implicit none
  complex(8),intent(in) :: zutmp(NL, NBoccmax,NK_s:NK_e)
  real(8),intent(in)    :: IaLxyz
  real(8),intent(out)   :: jx(NKB),jy(NKB),jz(NKB)
  integer,intent(in)    :: ikb_s,ikb_e

  integer    :: ilma,ia,j,i,ix,iy,iz
  real(8)    :: x,y,z
  real(8)    :: jxt,jyt,jzt
  complex(8) :: uVpsi,uVpsix,uVpsiy,uVpsiz
  integer    :: ikb, ik,ib

!$acc data pcopy(jx,jy,jz) pcopyin(zutmp) create(t4cp_uVpsix,t4cp_uVpsiy,t4cp_uVpsiz) &
!$acc& pcopyin(ik_table,ib_table,a_tbl,Mps,Jxyz,Jxx,Jyy,Jzz,lx,ly,lz,uV,ekr_omp,iuV,occ)
!$acc kernels
!$acc loop collapse(2) gang
  do ikb = ikb_s, ikb_e
    do ilma=1,Nlma
      ik=ik_table(ikb)
      ib=ib_table(ikb)
      ia=a_tbl(ilma)
      uVpsi=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
!$acc loop gang vector(256) reduction(+:uVpsi,uVpsix,uVpsiy,uVpsiz)
      do j=1,Mps(ia)
        i=Jxyz(j,ia)

        ix=Jxx(j,ia); x=Lx(i)*Hx-ix*aLx
        iy=Jyy(j,ia); y=Ly(i)*Hy-iy*aLy
        iz=Jzz(j,ia); z=Lz(i)*Hz-iz*aLz

        uVpsi =uVpsi +uV(j,ilma)*ekr_omp(j,ia,ik)  *zutmp(i,ib,ik)
        uVpsix=uVpsix+uV(j,ilma)*ekr_omp(j,ia,ik)*x*zutmp(i,ib,ik)
        uVpsiy=uVpsiy+uV(j,ilma)*ekr_omp(j,ia,ik)*y*zutmp(i,ib,ik)
        uVpsiz=uVpsiz+uV(j,ilma)*ekr_omp(j,ia,ik)*z*zutmp(i,ib,ik)
      end do
      uVpsi =uVpsi *Hxyz*iuV(ilma)
      uVpsix=uVpsix*Hxyz
      uVpsiy=uVpsiy*Hxyz
      uVpsiz=uVpsiz*Hxyz

      t4cp_uVpsix(ilma,ikb)=imag(conjg(uVpsix)*uVpsi)
      t4cp_uVpsiy(ilma,ikb)=imag(conjg(uVpsiy)*uVpsi)
      t4cp_uVpsiz(ilma,ikb)=imag(conjg(uVpsiz)*uVpsi)
    enddo
  enddo

!$acc loop gang
  do ikb = ikb_s, ikb_e
    ik=ik_table(ikb)
    ib=ib_table(ikb)

    jxt=0d0
    jyt=0d0
    jzt=0d0
!$acc loop vector(256) reduction(+:jxt,jyt,jzt)
    do ilma=1,Nlma
      jxt=jxt + t4cp_uVpsix(ilma,ikb)
      jyt=jyt + t4cp_uVpsiy(ilma,ikb)
      jzt=jzt + t4cp_uVpsiz(ilma,ikb)
    enddo
    jxt=jxt * occ(ib,ik)*IaLxyz*2
    jyt=jyt * occ(ib,ik)*IaLxyz*2
    jzt=jzt * occ(ib,ik)*IaLxyz*2

    jx(ikb)=jx(ikb)*Hxyz*IaLxyz + jxt
    jy(ikb)=jy(ikb)*Hxyz*IaLxyz + jyt
    jz(ikb)=jz(ikb)*Hxyz*IaLxyz + jzt
  enddo
!$acc end kernels
!$acc end data
end Subroutine
#endif ! ARTED_LBLK
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine current_GS_omp_KB_impl(zutmp,jxs,jys,jzs)
  use Global_Variables
  implicit none
  complex(8),intent(in) :: zutmp(0:NL-1,NB,NK_s:NK_e)
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
Subroutine current_omp_KB_impl(zutmp,jxs,jys,jzs)
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
#ifdef ARTED_CURRENT_OPTIMIZED
    call current_RT_stencil_impl(ik,ib,jx,jy,jz)
#else
    call current_stencil_impl(ik,ib,zutmp(:,ib,ik),nabt,jx,jy,jz)
#endif
    call current_pseudo_impl(ik,ib,zutmp(:,ib,ik),IaLxyz,jx,jy,jz)

    jxs=jxs+jx
    jys=jys+jy
    jzs=jzs+jz
  enddo
!$omp end parallel
end Subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_LBLK
subroutine current_acc_KB_impl(zutmp,jxs,jys,jzs)
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
#ifndef ARTED_CURRENT_OPTIMIZED
    call current_stencil_LBLK(zutmp(:,:,:), ikb_s,ikb_e)
#endif
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
Subroutine current_result_impl(jx,jy,jz)
  use Global_Variables
  use timelog
  implicit none
  real(8),intent(in) :: jx,jy,jz

  real(8) :: jav_l(3)

  call timelog_begin(LOG_ALLREDUCE)
  jav_l(1)=jx
  jav_l(2)=jy
  jav_l(3)=jz

  call MPI_ALLREDUCE(jav_l,jav,3,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)
  call timelog_end(LOG_ALLREDUCE)
end Subroutine
