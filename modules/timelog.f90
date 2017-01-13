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
module timelog
  implicit none

#ifdef ARTED_USE_TLOG
  include 'tlogf.h'

  integer,parameter :: LOG_DT_EVOLVE    = TLOG_EVENT_1_IN
  integer,parameter :: LOG_HPSI         = TLOG_EVENT_2_IN
  integer,parameter :: LOG_PSI_RHO      = TLOG_EVENT_3_IN
  integer,parameter :: LOG_HARTREE      = TLOG_EVENT_4_IN
  integer,parameter :: LOG_EXC_COR      = TLOG_EVENT_5_IN
  integer,parameter :: LOG_CURRENT      = TLOG_EVENT_6_IN
  integer,parameter :: LOG_TOTAL_ENERGY = TLOG_EVENT_7_IN
  integer,parameter :: LOG_ION_FORCE    = TLOG_EVENT_8_IN
  integer,parameter :: LOG_DT_EVOLVE_AC = TLOG_EVENT_9_IN
  integer,parameter :: LOG_K_SHIFT_WF   = TLOG_EVENT_10_IN
  integer,parameter :: LOG_ALLREDUCE    = TLOG_EVENT_11_IN
  integer,parameter :: LOG_OTHER        = TLOG_EVENT_12_IN

  integer,parameter :: LOG_CG           = TLOG_EVENT_13_IN
  integer,parameter :: LOG_DIAG         = TLOG_EVENT_14_IN
  integer,parameter :: LOG_SP_ENERGY    = TLOG_EVENT_15_IN
  integer,parameter :: LOG_GRAM_SCHMIDT = TLOG_EVENT_16_IN

  integer,parameter :: LOG_HPSI_INIT    = TLOG_EVENT_17_IN
  integer,parameter :: LOG_HPSI_STENCIL = TLOG_EVENT_18_IN
  integer,parameter :: LOG_HPSI_PSEUDO  = TLOG_EVENT_19_IN
  integer,parameter :: LOG_HPSI_UPDATE  = TLOG_EVENT_20_IN

  integer,parameter :: LOG_DYNAMICS     = TLOG_EVENT_21_IN

  integer,private,parameter :: LOG_OUT_OFFSET =TLOG_EVENT_1_OUT-TLOG_EVENT_1_IN
  integer,private,parameter :: LOG_ID_OFFSET  =TLOG_EVENT_1_IN
#else
  integer,parameter :: LOG_DT_EVOLVE    = 0
  integer,parameter :: LOG_HPSI         = 1
  integer,parameter :: LOG_PSI_RHO      = 2
  integer,parameter :: LOG_HARTREE      = 3
  integer,parameter :: LOG_EXC_COR      = 4
  integer,parameter :: LOG_CURRENT      = 5
  integer,parameter :: LOG_TOTAL_ENERGY = 6
  integer,parameter :: LOG_ION_FORCE    = 7
  integer,parameter :: LOG_DT_EVOLVE_AC = 8
  integer,parameter :: LOG_K_SHIFT_WF   = 9
  integer,parameter :: LOG_OTHER        = 10
  integer,parameter :: LOG_ALLREDUCE    = 11

  integer,parameter :: LOG_CG           = 12
  integer,parameter :: LOG_DIAG         = 13
  integer,parameter :: LOG_SP_ENERGY    = 14
  integer,parameter :: LOG_GRAM_SCHMIDT = 15

  integer,parameter :: LOG_HPSI_INIT    = 16
  integer,parameter :: LOG_HPSI_STENCIL = 17
  integer,parameter :: LOG_HPSI_PSEUDO  = 18
  integer,parameter :: LOG_HPSI_UPDATE  = 19

  integer,parameter :: LOG_DYNAMICS     = 20

  integer,private,parameter :: LOG_OUT_OFFSET =0
  integer,private,parameter :: LOG_ID_OFFSET  =0
#endif

  integer,private,parameter   :: LOG_SIZE = 30
  real(8),private,allocatable :: log_time(:)
  real(8),private,allocatable :: log_temp(:)

  real(8),private,allocatable :: log_time_t(:,:)
  real(8),private,allocatable :: log_temp_t(:,:)

#ifdef ARTED_USE_TLOG
  integer,private :: log_verbose = 0
#endif

contains
  subroutine timelog_initialize
    use omp_lib, only: omp_get_max_threads
    implicit none
    allocate(log_time(0:LOG_SIZE - 1))
    allocate(log_temp(0:LOG_SIZE - 1))
    allocate(log_time_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    allocate(log_temp_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    call timelog_reset
  end subroutine

  subroutine timelog_set(e,t)
    implicit none
    integer,intent(in) :: e
    real(8),intent(in) :: t
    log_time(e - LOG_ID_OFFSET) = t
    log_temp(e - LOG_ID_OFFSET) = 0.d0
  end subroutine

  subroutine timelog_reset(e)
    implicit none
    integer,intent(in),optional :: e
    if(present(e)) then
      log_time  (e - LOG_ID_OFFSET)   = 0.d0
      log_temp  (e - LOG_ID_OFFSET)   = 0.d0
      log_time_t(e - LOG_ID_OFFSET,:) = 0.d0
      log_temp_t(e - LOG_ID_OFFSET,:) = 0.d0
    else
      log_time  (:)   = 0.d0
      log_temp  (:)   = 0.d0
      log_time_t(:,:) = 0.d0
      log_time_t(:,:) = 0.d0
    end if
  end subroutine

  subroutine timelog_reentrance_read(fd)
    implicit none
    integer,intent(in) :: fd
    read(fd) log_time(0:LOG_SIZE - 1)
    read(fd) log_temp(0:LOG_SIZE - 1)
  end subroutine

  subroutine timelog_reentrance_write(fd)
    implicit none
    integer,intent(in) :: fd
    write(fd) log_time(0:LOG_SIZE - 1)
    write(fd) log_temp(0:LOG_SIZE - 1)
  end subroutine

  subroutine timelog_enable_verbose
    implicit none
#ifdef ARTED_USE_TLOG
    call tlog_initialize_
    log_verbose = 1
#endif
  end subroutine

  subroutine timelog_disable_verbose
    implicit none
#ifdef ARTED_USE_TLOG
    call tlog_finalize_
    log_verbose = 0
#endif
  end subroutine

  function get_wtime()
    implicit none
    real(8) :: get_wtime
    real(8) :: omp_get_wtime
    get_wtime = omp_get_wtime()
  end function

  subroutine timelog_begin(e)
    implicit none
    integer,intent(in) :: e
    integer :: id
#ifdef ARTED_USE_TLOG
    call tlog_log_(e)
#endif
    id = e - LOG_ID_OFFSET
    log_temp(id) = get_wtime()
  end subroutine

  subroutine timelog_end(e)
    implicit none
    integer,intent(in) :: e
    integer :: id
#ifdef ARTED_USE_TLOG
    call tlog_log_(e + LOG_OUT_OFFSET)
#endif
    id = e - LOG_ID_OFFSET
    log_time(id) = log_time(id) + get_wtime() - log_temp(id)
  end subroutine

  subroutine timelog_thread_begin(e)
    use omp_lib
    implicit none
    integer,intent(in) :: e
    integer :: tid,id
    tid = omp_get_thread_num()
    if (tid == 0) then
      call timelog_begin(e)
    end if
    id  = e - LOG_ID_OFFSET
    log_temp_t(id,tid) = get_wtime()
  end subroutine

  subroutine timelog_thread_end(e)
    use omp_lib
    implicit none
    integer,intent(in) :: e
    integer :: tid,id
    tid = omp_get_thread_num()
    if (tid == 0) then
      call timelog_end(e)
    end if
    id = e - LOG_ID_OFFSET
    log_time_t(id,tid) = log_time_t(id,tid) + get_wtime() - log_temp_t(id,tid)
  end subroutine

  subroutine timelog_show_hour(str, e)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: e
    real(8) :: time,hour
    time = log_time(e - LOG_ID_OFFSET)
    hour = time / 3600
    write(*,*) str,time,'sec =',hour,'hour'
  end subroutine

  subroutine timelog_show_min(str, e)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: e
    real(8) :: time,mini
    time = log_time(e - LOG_ID_OFFSET)
    mini = time / 60
    write(*,*) str,time,'sec =',mini,'min'
  end subroutine

  subroutine timelog_write(fd,str,e)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: fd,e
    real(8) :: time
    time = log_time(e - LOG_ID_OFFSET)
    write(fd,*) str,time,'sec'
  end subroutine

  subroutine timelog_thread_write(fd,str,e)
    use omp_lib
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: fd,e
    real(8) :: time
    integer :: i
    write(fd,*) str
    do i=0,omp_get_max_threads()-1
      time = log_time_t(e - LOG_ID_OFFSET,i)
      write(fd,*) 'tid =',i,': ',time,'sec'
    end do
  end subroutine

  function timelog_get(e)
    implicit none
    integer,intent(in) :: e
    real(8)            :: timelog_get
    timelog_get = log_time(e - LOG_ID_OFFSET)
  end function

  function timelog_thread_get(e,tid)
    implicit none
    integer,intent(in) :: e,tid
    real(8)            :: timelog_thread_get
    timelog_thread_get = log_time_t(e - LOG_ID_OFFSET,tid)
  end function
end module timelog
