!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of hpsi.f90.
!
!  hpsi.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  hpsi.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with hpsi.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!This file is "hpsi.f90"
!This file contain a subroutine.
!Subroutine hpsi(q)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine hpsi(ik)
  use Global_Variables
  implicit none
  integer :: i,ia,j,ilma,ik
  real(8) :: k2
  complex(8) :: uVpsi,zs(3),zt(3)
  real(8) :: k2lap0_2
  real(8) :: nabxt(4),nabyt(4),nabzt(4)

  k2=sum(kAc(ik,:)**2)



  k2lap0_2=(k2-(lapx(0)+lapy(0)+lapz(0)))*0.5d0
  nabxt(1:4)=kAc(ik,1)*nabx(1:4)
  nabyt(1:4)=kAc(ik,2)*naby(1:4)
  nabzt(1:4)=kAc(ik,3)*nabz(1:4)

  select case(Nd)
  case(4)
    do i=1,NL
      zs(1)=&!lapx(0)*tpsi(i)&
          &+lapx(1)*(tpsi(ifdx(1,i))+tpsi(ifdx(-1,i)))&
          &+lapx(2)*(tpsi(ifdx(2,i))+tpsi(ifdx(-2,i)))&
          &+lapx(3)*(tpsi(ifdx(3,i))+tpsi(ifdx(-3,i)))&
          &+lapx(4)*(tpsi(ifdx(4,i))+tpsi(ifdx(-4,i)))

      zt(1)=&!nabx(0)*tpsi(i)&
          &+nabxt(1)*(tpsi(ifdx(1,i))-tpsi(ifdx(-1,i)))&
          &+nabxt(2)*(tpsi(ifdx(2,i))-tpsi(ifdx(-2,i)))&
          &+nabxt(3)*(tpsi(ifdx(3,i))-tpsi(ifdx(-3,i)))&
          &+nabxt(4)*(tpsi(ifdx(4,i))-tpsi(ifdx(-4,i)))

      zs(2)=&!lapy(0)*tpsi(i)&
          &+lapy(1)*(tpsi(ifdy(1,i))+tpsi(ifdy(-1,i)))&
          &+lapy(2)*(tpsi(ifdy(2,i))+tpsi(ifdy(-2,i)))&
          &+lapy(3)*(tpsi(ifdy(3,i))+tpsi(ifdy(-3,i)))&
          &+lapy(4)*(tpsi(ifdy(4,i))+tpsi(ifdy(-4,i)))

      zt(2)=&!naby(0)*tpsi(i)&
          &+nabyt(1)*(tpsi(ifdy(1,i))-tpsi(ifdy(-1,i)))&
          &+nabyt(2)*(tpsi(ifdy(2,i))-tpsi(ifdy(-2,i)))&
          &+nabyt(3)*(tpsi(ifdy(3,i))-tpsi(ifdy(-3,i)))&
          &+nabyt(4)*(tpsi(ifdy(4,i))-tpsi(ifdy(-4,i)))

      zs(3)=&!lapz(0)*tpsi(i)&
          &+lapz(1)*(tpsi(ifdz(1,i))+tpsi(ifdz(-1,i)))&
          &+lapz(2)*(tpsi(ifdz(2,i))+tpsi(ifdz(-2,i)))&
          &+lapz(3)*(tpsi(ifdz(3,i))+tpsi(ifdz(-3,i)))&
          &+lapz(4)*(tpsi(ifdz(4,i))+tpsi(ifdz(-4,i)))

      zt(3)=&!nabz(0)*tpsi(i)&
          &+nabzt(1)*(tpsi(ifdz(1,i))-tpsi(ifdz(-1,i)))&
          &+nabzt(2)*(tpsi(ifdz(2,i))-tpsi(ifdz(-2,i)))&
          &+nabzt(3)*(tpsi(ifdz(3,i))-tpsi(ifdz(-3,i)))&
          &+nabzt(4)*(tpsi(ifdz(4,i))-tpsi(ifdz(-4,i)))

      ttpsi(i)=k2lap0_2*tpsi(i)-0.5d0*(zs(1)+zs(2)+zs(3))-zI*(zt(1)+zt(2)+zt(3))
    enddo
  case default
    call err_finalize('Nd /= 4')
  end select

  htpsi(1:NL)=ttpsi(1:NL)+Vloc(1:NL)*tpsi(1:NL)

!Calculating nonlocal part
  do ilma=1,Nlma
    ia=a_tbl(ilma)
    uVpsi=0.d0
    do j=1,Mps(ia)
      i=Jxyz(j,ia)
      uVpsi=uVpsi+uV(j,ilma)*ekr(j,ia)*tpsi(i)
    enddo
    uVpsi=uVpsi*Hxyz*iuV(ilma)
    do j=1,Mps(ia)
      i=Jxyz(j,ia)
      htpsi(i)=htpsi(i)+conjg(ekr(j,ia))*uVpsi*uV(j,ilma)
    enddo
  enddo

  return
End Subroutine hpsi
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine hpsi_omp_KB(ik,tpsi,ttpsi,htpsi)
  use Global_Variables, only: functional, NL
  implicit none
  integer,intent(in)     :: ik
  complex(8),intent(in)  :: tpsi(NL)
  complex(8),intent(out) :: ttpsi(NL),htpsi(NL)

  select case(functional)
    case('PZ','PBE','TBmBJ')
      call hpsi1(ik,tpsi,ttpsi,htpsi)
    case('TPSS','VS98')
      call err_finalize('hpsi_omp_KB: TPSS/VS98 ver. not implemented.')
  end select

contains
  subroutine hpsi1(ik,tpsi_,ttpsi_,htpsi_)
    use Global_Variables
    use opt_variables
    use timelog
    implicit none
    integer,intent(in)     :: ik
    complex(8),intent(in)  :: tpsi_(NL)
    complex(8),intent(out) :: ttpsi_(NL),htpsi_(NL)

    integer :: i,ia,j,ilma
    real(8) :: k2
    complex(8) :: uVpsi
    real(8) :: k2lap0_2
    real(8) :: nabt(12)

    call timelog_thread_begin(LOG_HPSI)

    k2=sum(kAc(ik,:)**2)
    k2lap0_2=(k2-(lapx(0)+lapy(0)+lapz(0)))*0.5d0

    call timelog_thread_begin(LOG_HPSI_STENCIL)
    nabt( 1: 4)=kAc(ik,1)*nabx(1:4)
    nabt( 5: 8)=kAc(ik,2)*naby(1:4)
    nabt( 9:12)=kAc(ik,3)*nabz(1:4)
    call hpsi1_tuned(k2lap0_2,Vloc,lapt,nabt,tpsi_,ttpsi_,htpsi_)
    call timelog_thread_end(LOG_HPSI_STENCIL)

    call timelog_thread_begin(LOG_HPSI_PSEUDO)

    !Calculating nonlocal part
    do ilma=1,Nlma
      ia=a_tbl(ilma)
      uVpsi=0.d0
      do j=1,Mps(ia)
        i=Jxyz(j,ia)
        uVpsi=uVpsi+uV(j,ilma)*ekr_omp(j,ia,ik)*tpsi_(i)
      enddo
      uVpsi=uVpsi*Hxyz*iuV(ilma)
!dir$ ivdep
      do j=1,Mps(ia)
        i=Jxyz(j,ia)
        htpsi_(i)=htpsi_(i)+conjg(ekr_omp(j,ia,ik))*uVpsi*uV(j,ilma)
      enddo
    enddo
    call timelog_thread_end(LOG_HPSI_PSEUDO)

    call timelog_thread_end(LOG_HPSI)
  end subroutine
end subroutine hpsi_omp_KB
!!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Subroutine hpsi_omp_KB(ik,thr_id)
!  use Global_Variables
!  implicit none
!  integer :: thr_id,ik
!  select case(functional)
!    case('PZ','PBE','TBmBJ')
!    call hpsi1(ik,thr_id)
!  case('TPSS','VS98')
!    call hpsi2(ik,thr_id)
!  end select
!end Subroutine hpsi_omp_KB
!!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Subroutine hpsi1(ik,thr_id)
!  use Global_Variables
!  use timelog
!  implicit none
!  integer :: thr_id
!  integer :: i,ia,j,ilma,ik
!  real(8) :: k2
!  complex(8) :: uVpsi,zs(3),zt(3)
!  real(8) :: k2lap0_2
!  real(8) :: nabxt(4),nabyt(4),nabzt(4)
!
!  k2=sum(kAc(ik,:)**2)
!  k2lap0_2=(k2-(lapx(0)+lapy(0)+lapz(0)))*0.5d0
!  nabxt(1:4)=kAc(ik,1)*nabx(1:4)
!  nabyt(1:4)=kAc(ik,2)*naby(1:4)
!  nabzt(1:4)=kAc(ik,3)*nabz(1:4)
!
!  select case(Nd)
!  case(4)
!    do i=1,NL
!      zs(1)=&!lapx(0)*tpsi(i)&
!          &+lapx(1)*(tpsi_omp(ifdx(1,i),thr_id)+tpsi_omp(ifdx(-1,i),thr_id))&
!          &+lapx(2)*(tpsi_omp(ifdx(2,i),thr_id)+tpsi_omp(ifdx(-2,i),thr_id))&
!          &+lapx(3)*(tpsi_omp(ifdx(3,i),thr_id)+tpsi_omp(ifdx(-3,i),thr_id))&
!          &+lapx(4)*(tpsi_omp(ifdx(4,i),thr_id)+tpsi_omp(ifdx(-4,i),thr_id))
!
!      zt(1)=&!nabx(0)*tpsi(i)&
!          &+nabxt(1)*(tpsi_omp(ifdx(1,i),thr_id)-tpsi_omp(ifdx(-1,i),thr_id))&
!          &+nabxt(2)*(tpsi_omp(ifdx(2,i),thr_id)-tpsi_omp(ifdx(-2,i),thr_id))&
!          &+nabxt(3)*(tpsi_omp(ifdx(3,i),thr_id)-tpsi_omp(ifdx(-3,i),thr_id))&
!          &+nabxt(4)*(tpsi_omp(ifdx(4,i),thr_id)-tpsi_omp(ifdx(-4,i),thr_id))
!
!      zs(2)=&!lapy(0)*tpsi(i)&
!          &+lapy(1)*(tpsi_omp(ifdy(1,i),thr_id)+tpsi_omp(ifdy(-1,i),thr_id))&
!          &+lapy(2)*(tpsi_omp(ifdy(2,i),thr_id)+tpsi_omp(ifdy(-2,i),thr_id))&
!          &+lapy(3)*(tpsi_omp(ifdy(3,i),thr_id)+tpsi_omp(ifdy(-3,i),thr_id))&
!          &+lapy(4)*(tpsi_omp(ifdy(4,i),thr_id)+tpsi_omp(ifdy(-4,i),thr_id))
!
!      zt(2)=&!naby(0)*tpsi(i)&
!          &+nabyt(1)*(tpsi_omp(ifdy(1,i),thr_id)-tpsi_omp(ifdy(-1,i),thr_id))&
!          &+nabyt(2)*(tpsi_omp(ifdy(2,i),thr_id)-tpsi_omp(ifdy(-2,i),thr_id))&
!          &+nabyt(3)*(tpsi_omp(ifdy(3,i),thr_id)-tpsi_omp(ifdy(-3,i),thr_id))&
!          &+nabyt(4)*(tpsi_omp(ifdy(4,i),thr_id)-tpsi_omp(ifdy(-4,i),thr_id))
!
!      zs(3)=&!lapz(0)*tpsi(i)&
!          &+lapz(1)*(tpsi_omp(ifdz(1,i),thr_id)+tpsi_omp(ifdz(-1,i),thr_id))&
!          &+lapz(2)*(tpsi_omp(ifdz(2,i),thr_id)+tpsi_omp(ifdz(-2,i),thr_id))&
!          &+lapz(3)*(tpsi_omp(ifdz(3,i),thr_id)+tpsi_omp(ifdz(-3,i),thr_id))&
!          &+lapz(4)*(tpsi_omp(ifdz(4,i),thr_id)+tpsi_omp(ifdz(-4,i),thr_id))
!
!      zt(3)=&!nabz(0)*tpsi(i)&
!          &+nabzt(1)*(tpsi_omp(ifdz(1,i),thr_id)-tpsi_omp(ifdz(-1,i),thr_id))&
!          &+nabzt(2)*(tpsi_omp(ifdz(2,i),thr_id)-tpsi_omp(ifdz(-2,i),thr_id))&
!          &+nabzt(3)*(tpsi_omp(ifdz(3,i),thr_id)-tpsi_omp(ifdz(-3,i),thr_id))&
!          &+nabzt(4)*(tpsi_omp(ifdz(4,i),thr_id)-tpsi_omp(ifdz(-4,i),thr_id))
!
!      ttpsi_omp(i,thr_id)=k2lap0_2*tpsi_omp(i,thr_id)-0.5d0*(zs(1)+zs(2)+zs(3))-zI*(zt(1)+zt(2)+zt(3))
!    enddo
!  case default
!    call err_finalize('Nd /= 4')
!  end select
!
!  htpsi_omp(1:NL,thr_id)=ttpsi_omp(1:NL,thr_id)+Vloc(1:NL)*tpsi_omp(1:NL,thr_id)
!
!!Calculating nonlocal part
!  do ilma=1,Nlma
!    ia=a_tbl(ilma)
!    uVpsi=0.d0
!    do j=1,Mps(ia)
!      i=Jxyz(j,ia)
!      uVpsi=uVpsi+uV(j,ilma)*ekr_omp(j,ia,ik)*tpsi_omp(i,thr_id)
!    enddo
!    uVpsi=uVpsi*Hxyz*iuV(ilma)
!    do j=1,Mps(ia)
!      i=Jxyz(j,ia)
!      htpsi_omp(i,thr_id)=htpsi_omp(i,thr_id)+conjg(ekr_omp(j,ia,ik))*uVpsi*uV(j,ilma)
!    enddo
!  enddo
!
!  return
!End Subroutine hpsi1
!!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!Subroutine hpsi2(ik,thr_id)
!  use Global_Variables
!  use timelog
!  implicit none
!  integer :: thr_id
!  integer :: i,ia,j,ilma,ik
!  real(8) :: k2
!  complex(8) :: uVpsi,zs(3),zt(3)
!  real(8) :: k2lap0_2
!  real(8) :: nabxt(4),nabyt(4),nabzt(4)
!  complex(8) :: gtpsi(NL,3)
!
!  k2=sum(kAc(ik,:)**2)
!  k2lap0_2=(k2-(lapx(0)+lapy(0)+lapz(0)))*0.5d0
!  nabxt(1:4)=kAc(ik,1)*nabx(1:4)
!  nabyt(1:4)=kAc(ik,2)*naby(1:4)
!  nabzt(1:4)=kAc(ik,3)*nabz(1:4)
!
!  select case(Nd)
!  case(4)
!    do i=1,NL
!      zs(1)=&!lapx(0)*tpsi(i)&
!          &+lapx(1)*(tpsi_omp(ifdx(1,i),thr_id)+tpsi_omp(ifdx(-1,i),thr_id))&
!          &+lapx(2)*(tpsi_omp(ifdx(2,i),thr_id)+tpsi_omp(ifdx(-2,i),thr_id))&
!          &+lapx(3)*(tpsi_omp(ifdx(3,i),thr_id)+tpsi_omp(ifdx(-3,i),thr_id))&
!          &+lapx(4)*(tpsi_omp(ifdx(4,i),thr_id)+tpsi_omp(ifdx(-4,i),thr_id))
!
!      zt(1)=&!nabx(0)*tpsi(i)&
!          &+nabx(1)*(tpsi_omp(ifdx(1,i),thr_id)-tpsi_omp(ifdx(-1,i),thr_id))&
!          &+nabx(2)*(tpsi_omp(ifdx(2,i),thr_id)-tpsi_omp(ifdx(-2,i),thr_id))&
!          &+nabx(3)*(tpsi_omp(ifdx(3,i),thr_id)-tpsi_omp(ifdx(-3,i),thr_id))&
!          &+nabx(4)*(tpsi_omp(ifdx(4,i),thr_id)-tpsi_omp(ifdx(-4,i),thr_id))
!
!      gtpsi(i,1)=zt(1)+zI*kAc(ik,1)*tpsi_omp(i,thr_id)
!      zt(1)=zt(1)*kAc(ik,1)
!
!      zs(2)=&!lapy(0)*tpsi(i)&
!          &+lapy(1)*(tpsi_omp(ifdy(1,i),thr_id)+tpsi_omp(ifdy(-1,i),thr_id))&
!          &+lapy(2)*(tpsi_omp(ifdy(2,i),thr_id)+tpsi_omp(ifdy(-2,i),thr_id))&
!          &+lapy(3)*(tpsi_omp(ifdy(3,i),thr_id)+tpsi_omp(ifdy(-3,i),thr_id))&
!          &+lapy(4)*(tpsi_omp(ifdy(4,i),thr_id)+tpsi_omp(ifdy(-4,i),thr_id))
!
!      zt(2)=&!naby(0)*tpsi(i)&
!          &+naby(1)*(tpsi_omp(ifdy(1,i),thr_id)-tpsi_omp(ifdy(-1,i),thr_id))&
!          &+naby(2)*(tpsi_omp(ifdy(2,i),thr_id)-tpsi_omp(ifdy(-2,i),thr_id))&
!          &+naby(3)*(tpsi_omp(ifdy(3,i),thr_id)-tpsi_omp(ifdy(-3,i),thr_id))&
!          &+naby(4)*(tpsi_omp(ifdy(4,i),thr_id)-tpsi_omp(ifdy(-4,i),thr_id))
!
!      gtpsi(i,2)=zt(2)+zI*kAc(ik,2)*tpsi_omp(i,thr_id)
!      zt(2)=zt(2)*kAc(ik,2)
!
!      zs(3)=&!lapz(0)*tpsi(i)&
!          &+lapz(1)*(tpsi_omp(ifdz(1,i),thr_id)+tpsi_omp(ifdz(-1,i),thr_id))&
!          &+lapz(2)*(tpsi_omp(ifdz(2,i),thr_id)+tpsi_omp(ifdz(-2,i),thr_id))&
!          &+lapz(3)*(tpsi_omp(ifdz(3,i),thr_id)+tpsi_omp(ifdz(-3,i),thr_id))&
!          &+lapz(4)*(tpsi_omp(ifdz(4,i),thr_id)+tpsi_omp(ifdz(-4,i),thr_id))
!
!      zt(3)=&!nabz(0)*tpsi(i)&
!          &+nabz(1)*(tpsi_omp(ifdz(1,i),thr_id)-tpsi_omp(ifdz(-1,i),thr_id))&
!          &+nabz(2)*(tpsi_omp(ifdz(2,i),thr_id)-tpsi_omp(ifdz(-2,i),thr_id))&
!          &+nabz(3)*(tpsi_omp(ifdz(3,i),thr_id)-tpsi_omp(ifdz(-3,i),thr_id))&
!          &+nabz(4)*(tpsi_omp(ifdz(4,i),thr_id)-tpsi_omp(ifdz(-4,i),thr_id))
!
!      gtpsi(i,3)=zt(3)+zI*kAc(ik,3)*tpsi_omp(i,thr_id)
!      zt(3)=zt(3)*kAc(ik,3)
!
!      ttpsi_omp(i,thr_id)=k2lap0_2*tpsi_omp(i,thr_id)-0.5d0*(zs(1)+zs(2)+zs(3))-zI*(zt(1)+zt(2)+zt(3))
!    enddo
!  case default
!    call err_finalize('Nd /= 4')
!  end select
!
!  htpsi_omp(1:NL,thr_id)=ttpsi_omp(1:NL,thr_id)+Vloc(1:NL)*tpsi_omp(1:NL,thr_id)
!
!  do i=1,NL
!    htpsi_omp(i,thr_id)=htpsi_omp(i,thr_id)+ &
! &    zI/2*(tjr(i,1)*gtpsi(i,1)+tjr(i,2)*gtpsi(i,2)+tjr(i,3)*gtpsi(i,3))
!  enddo
!  htpsi_omp(:,thr_id)=htpsi_omp(:,thr_id)+0.5d0*tjr2(:)*tpsi_omp(:,thr_id)
!
!  gtpsi(:,1)=tmass(:)*gtpsi(:,1)-zI*tjr(:,1)*tpsi_omp(i,thr_id)
!  gtpsi(:,2)=tmass(:)*gtpsi(:,2)-zI*tjr(:,2)*tpsi_omp(i,thr_id)
!  gtpsi(:,3)=tmass(:)*gtpsi(:,3)-zI*tjr(:,3)*tpsi_omp(i,thr_id)
!
!  do i=1,NL
!    htpsi_omp(i,thr_id)=htpsi_omp(i,thr_id)-0.5d0*( &
!          &+nabx(1)*(gtpsi(ifdx(1,i),1)-gtpsi(ifdx(-1,i),1))&
!          &+nabx(2)*(gtpsi(ifdx(2,i),1)-gtpsi(ifdx(-2,i),1))&
!          &+nabx(3)*(gtpsi(ifdx(3,i),1)-gtpsi(ifdx(-3,i),1))&
!          &+nabx(4)*(gtpsi(ifdx(4,i),1)-gtpsi(ifdx(-4,i),1))&
!          &+naby(1)*(gtpsi(ifdy(1,i),2)-gtpsi(ifdy(-1,i),2))&
!          &+naby(2)*(gtpsi(ifdy(2,i),2)-gtpsi(ifdy(-2,i),2))&
!          &+naby(3)*(gtpsi(ifdy(3,i),2)-gtpsi(ifdy(-3,i),2))&
!          &+naby(4)*(gtpsi(ifdy(4,i),2)-gtpsi(ifdy(-4,i),2))&
!          &+nabz(1)*(gtpsi(ifdz(1,i),3)-gtpsi(ifdz(-1,i),3))&
!          &+nabz(2)*(gtpsi(ifdz(2,i),3)-gtpsi(ifdz(-2,i),3))&
!          &+nabz(3)*(gtpsi(ifdz(3,i),3)-gtpsi(ifdz(-3,i),3))&
!          &+nabz(4)*(gtpsi(ifdz(4,i),3)-gtpsi(ifdz(-4,i),3))&
!    & +zI*(kAc(ik,1)*gtpsi(i,1)+kAc(ik,2)*gtpsi(i,2)+kAc(ik,3)*gtpsi(i,3) ))
!  enddo
!
!!Calculating nonlocal part
!  do ilma=1,Nlma
!    ia=a_tbl(ilma)
!    uVpsi=0.d0
!    do j=1,Mps(ia)
!      i=Jxyz(j,ia)
!      uVpsi=uVpsi+uV(j,ilma)*ekr_omp(j,ia,ik)*tpsi_omp(i,thr_id)
!    enddo
!    uVpsi=uVpsi*Hxyz*iuV(ilma)
!    do j=1,Mps(ia)
!      i=Jxyz(j,ia)
!      htpsi_omp(i,thr_id)=htpsi_omp(i,thr_id)+conjg(ekr_omp(j,ia,ik))*uVpsi*uV(j,ilma)
!    enddo
!  enddo
!
!  return
!End Subroutine hpsi2
!!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
