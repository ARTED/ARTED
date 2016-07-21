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
