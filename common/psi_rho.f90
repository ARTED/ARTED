!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of psi_rho.f90.
!
!  psi_rho.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  psi_rho.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with psi_rho.f90.  If not, see <http://www.gnu.org/licenses/>.
!
!This file is "psi_rho.f90"
!This file contain one subroutine.
!Subroutine psi_rho
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef __INTEL_COMPILER
# define OMP_SIMD simd
#else
# define OMP_SIMD
#endif

subroutine psi_rho_GS
  use global_variables, only: zu_GS,NB
  implicit none
  call psi_rho_impl(zu_GS,NB)
end subroutine

subroutine psi_rho_RT
  use global_variables, only: zu,NBoccmax
  implicit none
  call psi_rho_impl(zu,NBoccmax)
end subroutine

subroutine psi_rho_impl(zutmp,zu_NB)
  use global_variables
  use timelog
  use opt_variables
  implicit none
  integer,intent(in)    :: zu_NB
  complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)

  call timelog_begin(LOG_PSI_RHO)

  select case(Sym)
  case(1)
    call sym1(zutmp,zu_NB,rho_l)
  case(4)
    if(crystal_structure == 'diamond')then
      call sym4(zutmp,zu_NB,rho_l,rho_tmp1)
    else
      call err_finalize('Bad crystal structure')
    end if
  case(8)
    if(crystal_structure == 'diamond')then
      call sym8(zutmp,zu_NB,rho_l,rho_tmp1,rho_tmp2)
    else
      call err_finalize('Bad crystal structure')
    end if
  case default
    call err_finalize('Bad Symmetry')
  end select

  call timelog_begin(LOG_ALLREDUCE)
  call MPI_ALLREDUCE(rho_l,rho,NL,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)
  call timelog_end(LOG_ALLREDUCE)

  call timelog_end(LOG_PSI_RHO)

contains
  subroutine reduce(tid,zfac,zutmp,zu_NB)
    use global_variables
    use opt_variables, only: zrhotmp, roundup_pow2
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: tid
    real(8),intent(in)    :: zfac
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)

    integer :: ib,ik,i,mytid

    mytid = tid

#ifdef ARTED_REDUCE_FOR_MANYCORE
    zrhotmp(:,mytid)=0.d0

!$omp do private(ik,ib,i) collapse(2)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax
    do i=0,NL-1
      zrhotmp(i,mytid)=zrhotmp(i,mytid)+(zfac*occ(ib,ik))*abs(zutmp(i,ib,ik))**2
    end do
    end do
    end do
!$omp end do

    i = roundup_pow2(NUMBER_THREADS/2)
    do while(i > 0)
      if(mytid < i) then
        zrhotmp(0:NL-1,mytid) = zrhotmp(0:NL-1,mytid) + zrhotmp(0:NL-1,mytid + i)
      end if
      i = i/2
!$omp barrier
    end do
#else
    mytid = 0

!$omp single
    zrhotmp(:,mytid) = 0.d0
!$omp end single

    do ik=NK_s,NK_e
    do ib=1,NBoccmax
!$omp do private(i)
    do i=0,NL-1
      zrhotmp(i,mytid)=zrhotmp(i,mytid)+(zfac*occ(ib,ik))*abs(zutmp(i,ib,ik))**2
    end do
!$omp end do
    end do
    end do
#endif
  end subroutine

  subroutine sym1(zutmp,zu_NB,zrho_l)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    integer :: tid

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()
    call reduce(tid,1.0d0,zutmp,zu_NB)
!$omp end parallel

    zrho_l(:) = zrhotmp(0:NL-1,0)
  end subroutine

  !====== diamond(4) structure =========================!
  subroutine sym4(zutmp,zu_NB,zrho_l,zrhotmp1)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    real(8) :: zrhotmp1(0:NL-1)
    integer :: i,tid
    real(8) :: zfac

    zfac=1.0d0/4d0

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()
    call reduce(tid,zfac,zutmp,zu_NB)

! 1.T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp(i,0)+zrhotmp(itable_sym(1,i+1)-1,0)
    end do
!$omp end do OMP_SIMD

! 2.T_3*T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrho_l(i)=zrhotmp1(i)+zrhotmp1(itable_sym(2,i+1)-1)
    end do
!$omp end do OMP_SIMD
!$omp end parallel
  end subroutine

  !====== diamond(8) structure =========================!
  subroutine sym8(zutmp,zu_NB,zrho_l,zrhotmp1,zrhotmp2)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    real(8) :: zrhotmp1(0:NL-1)
    real(8) :: zrhotmp2(0:NL-1)
    integer :: i,tid
    real(8) :: zfac

    ! wk(ik)=8.0,(ikx==iky >. wk(ik)=4.0)
    zfac=1.0d0/32d0

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()
    call reduce(tid,zfac,zutmp,zu_NB)

! 1.T_4
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp(i,0)+zrhotmp(itable_sym(4,i+1)-1,0)
    end do
!$omp end do OMP_SIMD

! 2.T_3*T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp2(i)+zrhotmp2(itable_sym(5,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 2.T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp1(i)+zrhotmp1(itable_sym(3,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 2.T_1
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp2(i)+zrhotmp2(itable_sym(1,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 2.T_2
!$omp do OMP_SIMD
    do i=0,NL-1
      zrho_l(i)=zrhotmp1(i)+zrhotmp1(itable_sym(2,i+1)-1)
    end do
!$omp end do OMP_SIMD
!$omp end parallel
  end subroutine
end subroutine psi_rho_impl
