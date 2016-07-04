!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of hpsi_RT.f90.
!
!  hpsi_RT.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  hpsi_RT.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with hpsi_RT.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine hpsi_omp_KB_RT(ik,tpsi,htpsi)
  use Global_Variables, only: functional
  use opt_variables, only: PNL
  implicit none
  integer,intent(in)     :: ik
  complex(8),intent(in)  :: tpsi(0:PNL-1)
  complex(8),intent(out) :: htpsi(0:PNL-1)

  select case(functional)
    case('PZ','PBE','TBmBJ')
      call hpsi1(ik,tpsi,htpsi)
    case('TPSS','VS98')
      call err_finalize('hpsi_omp_KB_RT: TPSS/VS98 ver. not implemented.')
  end select

contains
  subroutine pseudo_pt(ik,tpsi,htpsi)
    use Global_Variables, only: Mps,uV,iuV,Hxyz,ekr_omp,Nlma,a_tbl
#ifdef ARTED_STENCIL_PADDING
    use opt_variables, only: zJxyz => zKxyz,PNL
#else
    use opt_variables, only: zJxyz,PNL
#endif
    implicit none
    integer,intent(in)     :: ik
    complex(8),intent(in)  :: tpsi(0:PNL-1)
    complex(8),intent(out) :: htpsi(0:PNL-1)

    integer    :: ilma,ia,j,i
    complex(8) :: uVpsi

    !Calculating nonlocal part
    do ilma=1,Nlma
      ia=a_tbl(ilma)
      uVpsi=0.d0
      do j=1,Mps(ia)
        i=zJxyz(j,ia)
        uVpsi=uVpsi+uV(j,ilma)*ekr_omp(j,ia,ik)*tpsi(i)
      enddo
      uVpsi=uVpsi*Hxyz*iuV(ilma)
!dir$ ivdep
      do j=1,Mps(ia)
        i=zJxyz(j,ia)
        htpsi(i)=htpsi(i)+conjg(ekr_omp(j,ia,ik))*uVpsi*uV(j,ilma)
      enddo
    enddo
  end subroutine

  subroutine hpsi1(ik,tpsi,htpsi)
    use Global_Variables, only: kAc,lapx,lapy,lapz,nabx,naby,nabz,Vloc
    use opt_variables, only: lapt,PNLx,PNLy,PNLz
#ifdef ARTED_SC
    use timelog
#endif
    implicit none
    integer,intent(in)     :: ik
    complex(8),intent(in)  :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8),intent(out) :: htpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)

    real(8)    :: k2,k2lap0_2
    real(8)    :: nabt(12)

#ifdef ARTED_SC
# define TIMELOG_BEG(id) call timelog_thread_begin(id)
# define TIMELOG_END(id) call timelog_thread_end(id)
#else
# define TIMELOG_BEG(id)
# define TIMELOG_END(id)
#endif

    k2=sum(kAc(ik,:)**2)
    k2lap0_2=(k2-(lapx(0)+lapy(0)+lapz(0)))*0.5d0

    TIMELOG_BEG(LOG_HPSI_STENCIL)
    nabt( 1: 4)=kAc(ik,1)*nabx(1:4)
    nabt( 5: 8)=kAc(ik,2)*naby(1:4)
    nabt( 9:12)=kAc(ik,3)*nabz(1:4)

#ifdef ARTED_STENCIL_ORIGIN
    call hpsi1_RT_stencil(k2lap0_2,Vloc,lapt(1:4),lapt(5:8),lapt(9:12),nabt(1:4),nabt(5:8),nabt(9:12),tpsi,htpsi)
#else
    call hpsi1_RT_stencil(k2lap0_2,Vloc,lapt,nabt,tpsi,htpsi)
#endif
    TIMELOG_END(LOG_HPSI_STENCIL)

    TIMELOG_BEG(LOG_HPSI_PSEUDO)
    call pseudo_pt(ik,tpsi,htpsi)
    TIMELOG_END(LOG_HPSI_PSEUDO)
  end subroutine
end subroutine hpsi_omp_KB_RT
