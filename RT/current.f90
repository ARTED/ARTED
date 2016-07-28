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
Subroutine current
  use Global_Variables
  implicit none
  integer :: ib,ik,i,j,ix,iy,iz,ia,ilma
  real(8) :: jx,jy,jz,jxt,jyt,jzt,x,y,z,kr,jav_l(3)
  complex(8) :: uVpsi,uVpsix,uVpsiy,uVpsiz,zs(3)

  jx=0.d0; jy=0.d0; jz=0.d0
  do ik=NK_s,NK_e
    do ib=1,NBoccmax
      jxt=0d0
      do i=1,NL
        jxt=jxt+real(zu(i,ib,ik))**2+aimag(zu(i,ib,ik))**2
      end do
      jyt=jxt ; jzt=jxt
      jx=jx+occ(ib,ik)*kAc(ik,1)*jxt
      jy=jy+occ(ib,ik)*kAc(ik,2)*jyt
      jz=jz+occ(ib,ik)*kAc(ik,3)*jzt
      
      jxt=0d0;jyt=0d0;jzt=0d0
      select case(Nd)
      case(4)
      do i=1,NL
        zs(1)=nabx(1)*(zu(ifdx(1,i),ib,ik)-zu(ifdx(-1,i),ib,ik))&
            &+nabx(2)*(zu(ifdx(2,i),ib,ik)-zu(ifdx(-2,i),ib,ik))&
            &+nabx(3)*(zu(ifdx(3,i),ib,ik)-zu(ifdx(-3,i),ib,ik))&
            &+nabx(4)*(zu(ifdx(4,i),ib,ik)-zu(ifdx(-4,i),ib,ik))
        zs(2)=naby(1)*(zu(ifdy(1,i),ib,ik)-zu(ifdy(-1,i),ib,ik))&
            &+naby(2)*(zu(ifdy(2,i),ib,ik)-zu(ifdy(-2,i),ib,ik))&
            &+naby(3)*(zu(ifdy(3,i),ib,ik)-zu(ifdy(-3,i),ib,ik))&
            &+naby(4)*(zu(ifdy(4,i),ib,ik)-zu(ifdy(-4,i),ib,ik))
        zs(3)=nabz(1)*(zu(ifdz(1,i),ib,ik)-zu(ifdz(-1,i),ib,ik))&
            &+nabz(2)*(zu(ifdz(2,i),ib,ik)-zu(ifdz(-2,i),ib,ik))&
            &+nabz(3)*(zu(ifdz(3,i),ib,ik)-zu(ifdz(-3,i),ib,ik))&
            &+nabz(4)*(zu(ifdz(4,i),ib,ik)-zu(ifdz(-4,i),ib,ik))
        jxt=jxt+imag(conjg(zu(i,ib,ik))*zs(1))
        jyt=jyt+imag(conjg(zu(i,ib,ik))*zs(2))
        jzt=jzt+imag(conjg(zu(i,ib,ik))*zs(3))
      enddo
      case default
        call err_finalize('Nd /= 4')
      end select
      jx=jx+jxt*occ(ib,ik)
      jy=jy+jyt*occ(ib,ik)
      jz=jz+jzt*occ(ib,ik)

    enddo
  enddo

  jx=jx*Hxyz/aLxyz
  jy=jy*Hxyz/aLxyz
  jz=jz*Hxyz/aLxyz

  jxt=0d0;jyt=0d0;jzt=0d0
  do ik=NK_s,NK_e
    do ia=1,NI
      do j=1,Mps(ia)
        i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
        kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
        ekr(j,ia)=exp(zI*kr)
      enddo
    end do
    do ib=1,NBoccmax
      do ilma=1,Nlma
        ia=a_tbl(ilma)
        uVpsi=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
        do j=1,Mps(ia)
          i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
          x=Lx(i)*Hx-ix*aLx
          y=Ly(i)*Hy-iy*aLy
          z=Lz(i)*Hz-iz*aLz
          uVpsi =uVpsi +uV(j,ilma)*ekr(j,ia)  *zu(i,ib,ik)
          uVpsix=uVpsix+uV(j,ilma)*ekr(j,ia)*x*zu(i,ib,ik)
          uVpsiy=uVpsiy+uV(j,ilma)*ekr(j,ia)*y*zu(i,ib,ik)
          uVpsiz=uVpsiz+uV(j,ilma)*ekr(j,ia)*z*zu(i,ib,ik)
        end do
        uVpsi =uVpsi *Hxyz*iuV(ilma)
        uVpsix=uVpsix*Hxyz
        uVpsiy=uVpsiy*Hxyz
        uVpsiz=uVpsiz*Hxyz
        jxt=jxt+occ(ib,ik)/aLxyz*2*imag(conjg(uVpsix)*uVpsi)
        jyt=jyt+occ(ib,ik)/aLxyz*2*imag(conjg(uVpsiy)*uVpsi)
        jzt=jzt+occ(ib,ik)/aLxyz*2*imag(conjg(uVpsiz)*uVpsi)
      enddo
    enddo
  enddo

  jx=jx+jxt
  jy=jy+jyt
  jz=jz+jzt

  jav_l(1)=jx
  jav_l(2)=jy
  jav_l(3)=jz

  call MPI_ALLREDUCE(jav_l,jav,3,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)

  return
end Subroutine current
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine current_GS
  use Global_Variables
  implicit none
  integer :: ib,ik,i,j,ix,iy,iz,ia,ilma
  real(8) :: jx(NL),jy(NL),jz(NL),x,y,z,kr,jav_l(3)
  complex(8) :: uVpsi,uVpsix,uVpsiy,uVpsiz,zs(3)

  jx=0.d0; jy=0.d0; jz=0.d0
  do ik=NK_s,NK_e
    do ib=1,NBoccmax
      do i=1,NL
        jx(i)=jx(i)+occ(ib,ik)*kAc(ik,1)*abs(zu_GS(i,ib,ik))**2
        jy(i)=jy(i)+occ(ib,ik)*kAc(ik,2)*abs(zu_GS(i,ib,ik))**2
        jz(i)=jz(i)+occ(ib,ik)*kAc(ik,3)*abs(zu_GS(i,ib,ik))**2
      enddo

      select case(Nd)
      case(4)
      do i=1,NL
        zs(1)=nabx(1)*(zu_GS(ifdx(1,i),ib,ik)-zu_GS(ifdx(-1,i),ib,ik))&
            &+nabx(2)*(zu_GS(ifdx(2,i),ib,ik)-zu_GS(ifdx(-2,i),ib,ik))&
            &+nabx(3)*(zu_GS(ifdx(3,i),ib,ik)-zu_GS(ifdx(-3,i),ib,ik))&
            &+nabx(4)*(zu_GS(ifdx(4,i),ib,ik)-zu_GS(ifdx(-4,i),ib,ik))
        zs(2)=naby(1)*(zu_GS(ifdy(1,i),ib,ik)-zu_GS(ifdy(-1,i),ib,ik))&
            &+naby(2)*(zu_GS(ifdy(2,i),ib,ik)-zu_GS(ifdy(-2,i),ib,ik))&
            &+naby(3)*(zu_GS(ifdy(3,i),ib,ik)-zu_GS(ifdy(-3,i),ib,ik))&
            &+naby(4)*(zu_GS(ifdy(4,i),ib,ik)-zu_GS(ifdy(-4,i),ib,ik))
        zs(3)=nabz(1)*(zu_GS(ifdz(1,i),ib,ik)-zu_GS(ifdz(-1,i),ib,ik))&
            &+nabz(2)*(zu_GS(ifdz(2,i),ib,ik)-zu_GS(ifdz(-2,i),ib,ik))&
            &+nabz(3)*(zu_GS(ifdz(3,i),ib,ik)-zu_GS(ifdz(-3,i),ib,ik))&
            &+nabz(4)*(zu_GS(ifdz(4,i),ib,ik)-zu_GS(ifdz(-4,i),ib,ik))
        jx(i)=jx(i)+occ(ib,ik)*imag(conjg(zu_GS(i,ib,ik))*zs(1))
        jy(i)=jy(i)+occ(ib,ik)*imag(conjg(zu_GS(i,ib,ik))*zs(2))
        jz(i)=jz(i)+occ(ib,ik)*imag(conjg(zu_GS(i,ib,ik))*zs(3))
      enddo
      case default
        call err_finalize('Nd /= 4')
      end select

    enddo
  enddo

  jav_l(1)=sum(jx)*Hxyz/aLxyz
  jav_l(2)=sum(jy)*Hxyz/aLxyz
  jav_l(3)=sum(jz)*Hxyz/aLxyz

  do ik=NK_s,NK_e
    do ia=1,NI
      do j=1,Mps(ia)
        i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
        kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
        ekr(j,ia)=exp(zI*kr)
      enddo
    enddo
    do ib=1,NBoccmax
      do ilma=1,Nlma
        ia=a_tbl(ilma)
        uVpsi=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
        do j=1,Mps(ia)
          i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
          x=Lx(i)*Hx-ix*aLx
          y=Ly(i)*Hy-iy*aLy
          z=Lz(i)*Hz-iz*aLz
          uVpsi =uVpsi +uV(j,ilma)*ekr(j,ia)  *zu_GS(i,ib,ik)
          uVpsix=uVpsix+uV(j,ilma)*ekr(j,ia)*x*zu_GS(i,ib,ik)
          uVpsiy=uVpsiy+uV(j,ilma)*ekr(j,ia)*y*zu_GS(i,ib,ik)
          uVpsiz=uVpsiz+uV(j,ilma)*ekr(j,ia)*z*zu_GS(i,ib,ik)
        end do
        uVpsi =uVpsi *Hxyz*iuV(ilma)
        uVpsix=uVpsix*Hxyz
        uVpsiy=uVpsiy*Hxyz
        uVpsiz=uVpsiz*Hxyz
        jav_l(1)=jav_l(1)+occ(ib,ik)/aLxyz*2*imag(conjg(uVpsix)*uVpsi)
        jav_l(2)=jav_l(2)+occ(ib,ik)/aLxyz*2*imag(conjg(uVpsiy)*uVpsi)
        jav_l(3)=jav_l(3)+occ(ib,ik)/aLxyz*2*imag(conjg(uVpsiz)*uVpsi)
      enddo
    enddo
  enddo
  call MPI_ALLREDUCE(jav_l,jav,3,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)

  return
End Subroutine Current_GS
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
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
Subroutine current_GS_omp_KB
  use Global_Variables
  use timelog
  implicit none
  real(8) :: jx,jy,jz

  call timelog_begin(LOG_CURRENT)
  call current_GS_omp_KB_impl(zu_GS,jx,jy,jz)
  call current_result_impl(jx,jy,jz)
  call timelog_end(LOG_CURRENT)
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

#ifdef ARTED_USE_NVTX
  call nvtxStartRange('current_omp_KB',2)
#endif
  call timelog_begin(LOG_CURRENT)
  call current_omp_KB_impl(zu,jx,jy,jz)
  call current_result_impl(jx,jy,jz)
  call timelog_end(LOG_CURRENT)
#ifdef ARTED_USE_NVTX
  call nvtxEndRange()
#endif
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
