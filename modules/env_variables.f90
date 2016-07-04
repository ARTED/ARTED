!
!  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
!  Copyright (C) 2016  ARTED developers
!
!  This file is part of env_variables.f90.
!
!  env_variables.f90 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  env_variables.f90 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with env_variables.f90.  If not, see <http://www.gnu.org/licenses/>.
!
module environment
  implicit none

  real(8) :: CPU_TASK_RATIO
  real(8) :: MIC_TASK_RATIO

  integer :: CPU_PROCESS_PER_NODE
  integer :: MIC_PROCESS_PER_NODE

  integer :: ENABLE_LOAD_BALANCER

  integer :: IS_SYMMETRIC_DEBUG

contains
  subroutine load_environments
    implicit none

    call get_cpu_task_ratio_internal(CPU_TASK_RATIO)
    MIC_TASK_RATIO = 2.0 - CPU_TASK_RATIO

    call get_cpu_ppn_internal(CPU_PROCESS_PER_NODE)
    call get_mic_ppn_internal(MIC_PROCESS_PER_NODE)

    call get_load_balancer_flag_internal(ENABLE_LOAD_BALANCER)

    call get_sym_debug_mode_internal(IS_SYMMETRIC_DEBUG)
  end subroutine
end module environment
