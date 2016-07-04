!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of hpsi_compiler_vec_original.f90.
!
!  hpsi_compiler_vec_original.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  hpsi_compiler_vec_original.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with hpsi_compiler_vec_original.f90.  If not, see <http://www.gnu.org/licenses/>.
!
subroutine hpsi1_RT_stencil(A,B,Cx,Cy,Cz,Dx,Dy,Dz,E,F)
  use global_variables, only: NL
  use opt_variables,    only: zifdx, zifdy, zifdz
  implicit none
  real(8),intent(in)     :: A
  real(8),intent(in)     :: B(0:NL-1)
  real(8),intent(in)     :: Cx(4),Cy(4),Cz(4)
  real(8),intent(in)     :: Dx(4),Dy(4),Dz(4)
  complex(8),intent(in)  :: E(0:NL-1)
  complex(8),intent(out) :: F(0:NL-1)
  complex(8),parameter   :: zI = (0.d0, 1.d0)

  integer    :: i
  complex(8) :: v(3), w(3)

#define IDX zifdx
#define IDY zifdy
#define IDZ zifdz

  do i=0,NL-1
    v(1)=Cx(1)*(E(IDX(1,i))+E(IDX(-1,i))) &
    &   +Cx(2)*(E(IDX(2,i))+E(IDX(-2,i))) &
    &   +Cx(3)*(E(IDX(3,i))+E(IDX(-3,i))) &
    &   +Cx(4)*(E(IDX(4,i))+E(IDX(-4,i)))

    w(1)=Dx(1)*(E(IDX(1,i))-E(IDX(-1,i))) &
    &   +Dx(2)*(E(IDX(2,i))-E(IDX(-2,i))) &
    &   +Dx(3)*(E(IDX(3,i))-E(IDX(-3,i))) &
    &   +Dx(4)*(E(IDX(4,i))-E(IDX(-4,i)))

    v(2)=Cy(1)*(E(IDY(1,i))+E(IDY(-1,i))) &
    &   +Cy(2)*(E(IDY(2,i))+E(IDY(-2,i))) &
    &   +Cy(3)*(E(IDY(3,i))+E(IDY(-3,i))) &
    &   +Cy(4)*(E(IDY(4,i))+E(IDY(-4,i)))

    w(2)=Dy(1)*(E(IDY(1,i))-E(IDY(-1,i))) &
    &   +Dy(2)*(E(IDY(2,i))-E(IDY(-2,i))) &
    &   +Dy(3)*(E(IDY(3,i))-E(IDY(-3,i))) &
    &   +Dy(4)*(E(IDY(4,i))-E(IDY(-4,i)))

    v(3)=Cz(1)*(E(IDZ(1,i))+E(IDZ(-1,i))) &
    &   +Cz(2)*(E(IDZ(2,i))+E(IDZ(-2,i))) &
    &   +Cz(3)*(E(IDZ(3,i))+E(IDZ(-3,i))) &
    &   +Cz(4)*(E(IDZ(4,i))+E(IDZ(-4,i)))

    w(3)=Dz(1)*(E(IDZ(1,i))-E(IDZ(-1,i))) &
    &   +Dz(2)*(E(IDZ(2,i))-E(IDZ(-2,i))) &
    &   +Dz(3)*(E(IDZ(3,i))-E(IDZ(-3,i))) &
    &   +Dz(4)*(E(IDZ(4,i))-E(IDZ(-4,i)))

    F(i) = B(i)*E(i) + A*E(i) - 0.5d0*(v(1)+v(2)+v(3)) - zI*(w(1)+w(2)+w(3))
  end do
end subroutine
