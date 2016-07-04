!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of hpsi_compiler_vec_simple.f90.
!
!  hpsi_compiler_vec_simple.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  hpsi_compiler_vec_simple.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with hpsi_compiler_vec_simple.f90.  If not, see <http://www.gnu.org/licenses/>.
!
subroutine hpsi1_RT_stencil(A,B,C,D,E,F)
  use global_variables, only: NL
  use opt_variables,    only: zifdt
  implicit none
  real(8),intent(in)     :: A
  real(8),intent(in)     :: B(0:NL-1)
  real(8),intent(in)     :: C(12)
  real(8),intent(in)     :: D(12)
  complex(8),intent(in)  :: E(0:NL-1)
  complex(8),intent(out) :: F(0:NL-1)
  complex(8),parameter   :: zI = (0.d0, 1.d0)

  integer    :: i,j
  complex(8) :: v,w

#define IDX zifdt

#ifdef __INTEL_COMPILER
!dir$ vector nontemporal(F)
#endif
#ifdef __FUJITSU
!OCL simd
!OCL noalias
#endif
  do i=0,NL-1
    v = 0
    w = 0

    do j=1,12
      v = v + C(j) * (E(IDX(j,i)) + E(IDX(j+12,i)))
      w = w + D(j) * (E(IDX(j,i)) - E(IDX(j+12,i)))
    end do

    F(i) = B(i)*E(i) + A*E(i) - 0.5d0*v - zI*w
  end do
end subroutine
