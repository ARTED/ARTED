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
module performance_analyzer
  implicit none

  public  print_stencil_size, show_performance

  private summation_threads, get_gflops, get_hamiltonian_chunk_size
  private get_stencil_FLOP, get_pseudo_pt_FLOP, get_update_FLOP

contains
  subroutine print_stencil_size
    use global_variables, only: NK_s,NK_e,NBoccmax,NL,Nt,NUMBER_THREADS
    implicit none
    integer :: NK, NB

    NK = NK_e - NK_s + 1
    NB = NBoccmax
    print *, 'NK =', NK
    print *, 'NB =', NB
    print *, 'NL =', NL
    print *, 'Nt =', (Nt + 1)
    print *, 'Number of Domain/Thread =', real(NK * NB) / NUMBER_THREADS
  end subroutine

  subroutine show_performance
    use global_variables
    use timelog
    implicit none
    real(8) :: lgflops(4), tgflops(4)
#ifdef ARTED_MS
    real(8) :: sgflops(4)
#endif

    call summation_threads(lgflops)

#ifdef ARTED_MS
    call MPI_ALLREDUCE(lgflops,sgflops,4,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)
#endif
    call MPI_ALLREDUCE(lgflops,tgflops,4,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(Myrank == 0) then
      print *, '###'
      print *, 'Performance [GFLOPS]'
      print *, '  Processor'
      print *, '    hamiltonian  :', lgflops(4)
      print *, '    - stencil    :', lgflops(1)
      print *, '    - pseudo pt. :', lgflops(2)
      print *, '    - update     :', lgflops(3)
#ifdef ARTED_MS
      print *, '  Macroscopic grid point'
      print *, '    hamiltonian  :', sgflops(4)
      print *, '    - stencil    :', sgflops(1)
      print *, '    - pseudo pt. :', sgflops(2)
      print *, '    - update     :', sgflops(3)
#endif
      print *, '  System total'
      print *, '    hamiltonian  :', tgflops(4)
      print *, '    - stencil    :', tgflops(1)
      print *, '    - pseudo pt. :', tgflops(2)
      print *, '    - update     :', tgflops(3)
      print *, '###'
    endif
  end subroutine

  subroutine summation_threads(lgflops)
    use global_variables, only: NUMBER_THREADS
    use timelog
    implicit none
    real(8), intent(out) :: lgflops(4)
    real(8) :: hflop(3), htime(4)
    integer :: i, cnt
    integer :: chunk_size(0:NUMBER_THREADS-1)

    call get_hamiltonian_chunk_size(chunk_size)

    lgflops = 0.0d0
    do i=0,NUMBER_THREADS-1
      cnt = chunk_size(i)

      hflop(1) = get_stencil_FLOP(cnt)
      hflop(2) = get_pseudo_pt_FLOP(cnt)
      hflop(3) = get_update_FLOP(cnt)

      htime(1) = timelog_thread_get(LOG_HPSI_STENCIL, i)
      htime(2) = timelog_thread_get(LOG_HPSI_PSEUDO, i)
      htime(3) = timelog_thread_get(LOG_HPSI_UPDATE, i)
      htime(4) = timelog_thread_get(LOG_HPSI_INIT, i)

      lgflops(1) = lgflops(1) + get_gflops(hflop(1), htime(1))
      lgflops(2) = lgflops(2) + get_gflops(hflop(2), htime(2))
      lgflops(3) = lgflops(3) + get_gflops(hflop(3), htime(3))
      lgflops(4) = lgflops(4) + get_gflops(sum(hflop), sum(htime))
    end do
  end subroutine

  subroutine get_hamiltonian_chunk_size(chunk_size)
    use global_variables, only: NKB,NUMBER_THREADS
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(out) :: chunk_size(0:NUMBER_THREADS-1)
    integer :: ikb, tid

    tid = 0
    chunk_size(:) = 0
!$omp parallel private(tid) shared(chunk_size)
!$ tid = omp_get_thread_num()
!$omp do private(ikb)
    do ikb=1,NKB
      chunk_size(tid) = chunk_size(tid) + 1
    end do
!$omp end do
!$omp end parallel
  end subroutine

  function get_stencil_FLOP(chunk_size)
    use global_variables, only: NK_s,NK_e,NBoccmax,NL,Nt
    integer,intent(in),optional :: chunk_size
    real(8),parameter           :: FLOP = 158

    real(8) :: get_stencil_FLOP
    integer :: nsize

    if(present(chunk_size)) then
      nsize = chunk_size
    else
      nsize = (NK_e - NK_s + 1) * NBoccmax
    end if
    get_stencil_FLOP = nsize * 4*FLOP*NL * (Nt + 1)
  end function

  function get_pseudo_pt_FLOP(chunk_size)
    use global_variables, only: NK_s,NK_e,NBoccmax,Nt,a_tbl,Mps
    implicit none
    integer,intent(in),optional :: chunk_size
    real(8),parameter           :: FLOP_reduction = (2 + 6)     + 2
    real(8),parameter           :: FLOP_scatter   = (2 + 6 + 1) + 2 ! 1 = conjg(z)
    real(8),parameter           :: FLOP_scalar    = 2 + 2

    real(8) :: get_pseudo_pt_FLOP
    real(8) :: FLOP
    integer :: nsize

    FLOP = FLOP_scalar + (FLOP_reduction + FLOP_scatter) * sum(Mps(a_tbl(:)))

    if(present(chunk_size)) then
      nsize = chunk_size
    else
      nsize = (NK_e - NK_s + 1) * NBoccmax
    endif
    get_pseudo_pt_FLOP = nsize * 4*FLOP * (Nt + 1)
  end function

  function get_update_FLOP(chunk_size)
    use global_variables, only: NK_s,NK_e,NBoccmax,NL,Nt
    implicit none
    integer,intent(in),optional :: chunk_size
    real(8),parameter           :: FLOP = 6 + 2

    real(8) :: get_update_FLOP
    integer :: nsize

    if(present(chunk_size)) then
      nsize = chunk_size
    else
      nsize = (NK_e - NK_s + 1) * NBoccmax
    endif
    get_update_FLOP = nsize * 4*FLOP*NL * (Nt + 1)
  end function

  function get_gflops(FLOP,time)
    implicit none
    real(8),intent(in) :: FLOP
    real(8),intent(in) :: time
    real(8) :: get_gflops
    get_gflops = FLOP / (time * (10**9))
  end function
end module
