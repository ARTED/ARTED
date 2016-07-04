!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of current_compiler_vec.f90.
!
!  current_compiler_vec.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  current_compiler_vec.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with current_compiler_vec.f90.  If not, see <http://www.gnu.org/licenses/>.
!
subroutine current_stencil(C,E,F,G,H)
  use global_variables, only: NLx,NLy,NLz,zI
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  use opt_variables, only: LBX => STENCIL_BLOCKING_X, LBY => STENCIL_BLOCKING_Y
#endif
#ifndef ARTED_DOMAIN_POWER_OF_TWO
  use opt_variables, only: modx, mody, modz
#endif
  implicit none
  real(8),    intent(in)  :: C(12)
  complex(8), intent(in)  :: E(0:NLz-1,0:NLy-1,0:NLx-1)
  real(8),    intent(out) :: F,G,H

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  integer    :: bx,by
#endif
  integer    :: ix,iy,iz
  complex(8) :: v,w

#ifdef __INTEL_COMPILER
#if defined(__KNC__) || defined(__AVX512F__)
#   define MEM_ALIGN   64
#   define VECTOR_SIZE 4
# else
#   define MEM_ALIGN   32
#   define VECTOR_SIZE 2
# endif

!dir$ assume_aligned E:MEM_ALIGN
#endif

#ifdef ARTED_DOMAIN_POWER_OF_TWO
# ifdef __INTEL_COMPILER
!dir$ assume (mod(NLx, VECTOR_SIZE) == 0)
!dir$ assume (mod(NLy, VECTOR_SIZE) == 0)
!dir$ assume (mod(NLz, VECTOR_SIZE) == 0)
# endif
# define IDX(dt) iz,iy,and(ix+(dt)+NLx,NLx-1)
# define IDY(dt) iz,and(iy+(dt)+NLy,NLy-1),ix
# define IDZ(dt) and(iz+(dt)+NLz,NLz-1),iy,ix
#else
# define IDX(dt) iz,iy,modx(ix+(dt)+NLx)
# define IDY(dt) iz,mody(iy+(dt)+NLy),ix
# define IDZ(dt) modz(iz+(dt)+NLz),iy,ix
#endif

  H = 0
  G = 0
  F = 0

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  do bx=0,NLx-1,LBX
  do by=0,NLy-1,LBY
  do ix=bx,min(bx+LBX-1,NLx-1)
  do iy=by,min(by+LBY-1,NLy-1)
#else
  do ix=0,NLx-1
  do iy=0,NLy-1
#endif
#ifdef __INTEL_COMPILER
!dir$ simd
#endif
#ifdef __FUJITSU
!OCL simd
!OCL noalias
#endif
  do iz=0,NLz-1
    w = conjg(E(iz,iy,ix))

    v=(C( 9)*(E(IDZ(1))-E(IDZ(-1))) &
    & +C(10)*(E(IDZ(2))-E(IDZ(-2))) &
    & +C(11)*(E(IDZ(3))-E(IDZ(-3))) &
    & +C(12)*(E(IDZ(4))-E(IDZ(-4))))

    H = H + imag(w * v)

    v=(C( 5)*(E(IDY(1))-E(IDY(-1))) &
    & +C( 6)*(E(IDY(2))-E(IDY(-2))) &
    & +C( 7)*(E(IDY(3))-E(IDY(-3))) &
    & +C( 8)*(E(IDY(4))-E(IDY(-4))))

    G = G + imag(w * v)

    v=(C( 1)*(E(IDX(1))-E(IDX(-1))) &
    & +C( 2)*(E(IDX(2))-E(IDX(-2))) &
    & +C( 3)*(E(IDX(3))-E(IDX(-3))) &
    & +C( 4)*(E(IDX(4))-E(IDX(-4))))

    F = F + imag(w * v)
  end do
  end do
  end do
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  end do
  end do
#endif
end subroutine
