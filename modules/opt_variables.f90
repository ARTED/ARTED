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
module opt_variables
  implicit none

  real(8) :: lapt(12)

  integer                :: PNLx,PNLy,PNLz,PNL
  complex(8),allocatable :: ztpsi(:,:,:)

  real(8),allocatable :: zrhotmp(:,:)

  integer,allocatable :: zJxyz(:,:),zKxyz(:,:)

  integer,allocatable :: modx(:),mody(:),modz(:)

  integer :: STENCIL_BLOCKING_X
  integer :: STENCIL_BLOCKING_Y

  real(8),allocatable :: zcx(:,:),zcy(:,:),zcz(:,:)

#ifdef ARTED_STENCIL_ORIGIN
  integer,allocatable :: zifdx(:,:),zifdy(:,:),zifdz(:,:)
#endif

#ifndef ARTED_STENCIL_OPTIMIZED
  integer,allocatable :: zifdt(:,:)
#endif

  integer,allocatable :: hpsi_called(:)

#if defined(__KNC__) || defined(__AVX512F__)
# define MEM_ALIGNED 64
#else
# define MEM_ALIGNED 32
#endif

!dir$ attributes align:MEM_ALIGNED :: lapt
!dir$ attributes align:MEM_ALIGNED :: ztpsi
!dir$ attributes align:MEM_ALIGNED :: zrhotmp
!dir$ attributes align:MEM_ALIGNED :: zJxyz,zKxyz
!dir$ attributes align:MEM_ALIGNED :: zcx,zcy,zcz
!dir$ attributes align:MEM_ALIGNED :: modx,mody,modz

#ifdef ARTED_STENCIL_ORIGIN
!dir$ attributes align:MEM_ALIGNED :: zifdx,zifdy,zifdz
#endif

#ifndef ARTED_STENCIL_OPTIMIZED
!dir$ attributes align:MEM_ALIGNED :: zifdt
#endif

contains
  function ceil_power_of_two(n)
    implicit none
    integer,intent(in) :: n
    integer            :: ceil_power_of_two
    integer :: x
    x = n
    x = ior(x, ishft(x, 1))
    x = ior(x, ishft(x, 2))
    x = ior(x, ishft(x, 4))
    x = ior(x, ishft(x, 8))
    x = ior(x, ishft(x, 16))
    ceil_power_of_two = x - ishft(x, 1)
  end function

  function roundup_pow2(n)
    implicit none
    integer,intent(in) :: n
    integer            :: roundup_pow2,k

    k = n - 1
    k = ior(k, ishft(k,-1))
    k = ior(k, ishft(k,-2))
    k = ior(k, ishft(k,-4))
    k = ior(k, ishft(k,-8))
    k = ior(k, ishft(k,-16))

    roundup_pow2 = k + 1
  end function roundup_pow2

  subroutine opt_vars_initialize_p1
    use global_variables
    implicit none
    integer :: tid_range

    select case(functional)
      case('TPSS','VS98')
        call err_finalize('functional: TPSS/VS98 versions not implemented.')
    end select

#ifdef ARTED_REDUCE_FOR_MANYCORE
    tid_range = roundup_pow2(NUMBER_THREADS) - 1
#else
    tid_range = 0
#endif
    allocate(zrhotmp(0:NL-1,0:tid_range))
  end subroutine

  subroutine opt_vars_initialize_p2
    use global_variables
    implicit none

    integer :: ix,iy,iz

    PNLx = NLx
#ifdef ARTED_STENCIL_PADDING
    PNLy = NLy + 1
#else
    PNLy = NLy
#endif
    PNLz = NLz
    PNL  = PNLx * PNLy * PNLz

    allocate(ztpsi(0:PNL-1,4,0:NUMBER_THREADS-1))

    allocate(zcx(NBoccmax,NK_s:NK_e))
    allocate(zcy(NBoccmax,NK_s:NK_e))
    allocate(zcz(NBoccmax,NK_s:NK_e))

    lapt( 1: 4)=lapx(1:4)
    lapt( 5: 8)=lapy(1:4)
    lapt( 9:12)=lapz(1:4)

#ifdef ARTED_STENCIL_ORIGIN
    allocate(zifdx(-4:4,0:NL-1))
    allocate(zifdy(-4:4,0:NL-1))
    allocate(zifdz(-4:4,0:NL-1))

    zifdx(-4:4,0:NL-1) = ifdx(-4:4,1:NL) - 1
    zifdy(-4:4,0:NL-1) = ifdy(-4:4,1:NL) - 1
    zifdz(-4:4,0:NL-1) = ifdz(-4:4,1:NL) - 1
#endif

#ifndef ARTED_STENCIL_OPTIMIZED
    allocate(zifdt(24,0:NL-1))

    zifdt( 1: 4, 0:NL-1) = ifdx( 1: 4, 1:NL) - 1
    zifdt( 5: 8, 0:NL-1) = ifdy( 1: 4, 1:NL) - 1
    zifdt( 9:12, 0:NL-1) = ifdz( 1: 4, 1:NL) - 1
    zifdt(13:16, 0:NL-1) = ifdx(-1:-4:-1, 1:NL) - 1
    zifdt(17:20, 0:NL-1) = ifdy(-1:-4:-1, 1:NL) - 1
    zifdt(21:24, 0:NL-1) = ifdz(-1:-4:-1, 1:NL) - 1
#endif

    allocate(zJxyz(Nps,NI))

    zJxyz(1:Nps,1:NI) = Jxyz(1:Nps,1:NI) - 1

    allocate(modx(0:NLx*2+Nd-1))
    allocate(mody(0:NLy*2+Nd-1))
    allocate(modz(0:NLz*2+Nd-1))

    do ix=0,NLx*2+Nd-1
      modx(ix) = mod(ix,NLx)
    end do
    do iy=0,NLy*2+Nd-1
      mody(iy) = mod(iy,NLy)
    end do
    do iz=0,NLz*2+Nd-1
      modz(iz) = mod(iz,NLz)
    end do

    allocate(hpsi_called(0:NUMBER_THREADS-1))
    hpsi_called(:) = 0

#ifdef ARTED_STENCIL_PADDING
    call init_for_padding
#endif

#ifdef ARTED_STENCIL_LOOP_BLOCKING
    call auto_blocking
#endif
  end subroutine

  subroutine init_for_padding
    use global_variables
    implicit none
    integer :: a,ik,ix,iy,iz,jx,jy,jz,i,j,k
    real(8) :: x,y,z,r

    allocate(zKxyz(Nps,NI))

    do a=1,NI
      ik=Kion(a)
      j=0
      do ix=-2,2
      do iy=-2,2
      do iz=-2,2
        do jx=0,NLx-1
        do jy=0,NLy-1
        do jz=0,NLz-1
          i=jx* NLy* NLz + jy* NLz + jz + 1
          k=jx*PNLy*PNLz + jy*PNLz + jz
          x=Lx(i)*Hx-(Rion(1,a)+ix*aLx)
          y=Ly(i)*Hy-(Rion(2,a)+iy*aLy)
          z=Lz(i)*Hz-(Rion(3,a)+iz*aLz)
          r=sqrt(x*x+y*y+z*z)
          if (r<Rps(ik)) then
            j=j+1
            if(j<=Nps) then
              zKxyz(j,a)=k
            end if
          end if
        end do
        end do
        end do
      end do
      end do
      end do
    end do
  end subroutine

  subroutine auto_blocking
    implicit none
    integer,parameter :: L1cache_size =  8 * 1024
    integer,parameter :: value_size   = 24
    real(8) :: nyx
    integer :: sq

    nyx = dble(L1cache_size) / (PNLz * value_size)
    sq  = int(floor(sqrt(nyx)))

    STENCIL_BLOCKING_X = ceil_power_of_two(min(sq, PNLx))
    STENCIL_BLOCKING_Y = ceil_power_of_two(min(sq, PNLy))
  end subroutine

  subroutine symmetric_load_balancing(NK,NK_ave,NK_s,NK_e,NK_remainder,Myrank,Nprocs)
    use environment
    implicit none
    integer,intent(in)    :: NK
    integer,intent(in)    :: NK_ave
    integer,intent(inout) :: NK_s
    integer,intent(inout) :: NK_e
    integer,intent(inout) :: NK_remainder
    integer,intent(in)    :: Myrank
    integer,intent(in)    :: Nprocs

    integer :: NScpu,NSmic,NPcpu,NPmic,NPtotal
    integer :: np,npr,pos

    NPcpu   = CPU_PROCESS_PER_NODE
    NPmic   = MIC_PROCESS_PER_NODE
    NPtotal = NPcpu + NPmic

    if (Myrank == 0 .and. CPU_PROCESS_PER_NODE /= MIC_PROCESS_PER_NODE) then
      call err_finalize('CPU_PROCESS_PER_NODE /= MIC_PROCESS_PER_NODE')
    end if

    ! NK
    NScpu = int(NK_ave * CPU_TASK_RATIO)
    NSmic = int(NK_ave * MIC_TASK_RATIO)

    NK_remainder = NK - (NScpu * (Nprocs/2) + NSmic * (Nprocs/2))

    np  = myrank / NPtotal * NPtotal
    npr = mod(myrank, NPtotal)
    pos = (np / 2) * NScpu &
    &   + (np / 2) * NSmic

    if (npr < NPcpu) then
      pos = pos + npr * NScpu
    else
      pos = pos + NPcpu          * NScpu
      pos = pos + mod(npr,NPmic) * NSmic
    end if
    NK_s = pos
#ifdef __MIC__
    NK_e = pos + NSmic
#else
    NK_e = pos + NScpu
#endif
    if (Myrank+1 == Nprocs .and. NK_remainder /= 0) then
      NK_e = NK_e + NK_remainder
    end if
    NK_s = NK_s + 1

    ! Error check
    if(Myrank == Nprocs-1 .and. NK_e /= NK) then
      call err_finalize('prep. NK_e error')
    end if
  end subroutine

  function is_symmetric_mode()
    use global_variables
    implicit none
    integer             :: is_symmetric_mode
    integer             :: arch, ret, i
    integer,allocatable :: results(:)

    allocate(results(Nprocs))

#ifdef __MIC__
    arch = 1
#else
    arch = 2
#endif

    call MPI_Allgather(arch,1,MPI_INT,results,1,MPI_INT,MPI_COMM_WORLD,ierr)

    do i=2,Nprocs
      if(results(1) /= results(i)) then
        ret = 1
        exit
      end if
    end do

    is_symmetric_mode = ret
  end function

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

  function get_stencil_gflops(time,chunk_size)
    use global_variables, only: NK_s,NK_e,NBoccmax,NL,Nt
    implicit none
    real(8),intent(in)          :: time
    integer,intent(in),optional :: chunk_size
    real(8),parameter           :: FLOPS = 158

    real(8) :: get_stencil_gflops
    integer :: NK, nsize

    NK = NK_e - NK_s + 1
    if(present(chunk_size)) then
      nsize = chunk_size*NL
    else
      nsize = NK*NBoccmax*NL
    end if
    get_stencil_gflops = (nsize * 4*FLOPS * (Nt + 1)) / (time * (10**9))
  end function

  subroutine write_threads_performance
    use global_variables, only: directory,SYSname,NUMBER_THREADS
    use timelog
    implicit none
    integer,parameter :: fd = 32
    character(100)    :: filepath
    real(8)           :: time, lgflops, gflops
    integer           :: i, cnt
    filepath = trim(directory)//trim(SYSname)//'_hpsi_perf.log'

    open(fd,file=filepath,status='replace')
    call timelog_write(fd, '### global hpsi time', LOG_HPSI)
    call timelog_thread_write(fd, '### init time', LOG_HPSI_INIT)
    call timelog_thread_write(fd, '### stencil time', LOG_HPSI_STENCIL)
    call timelog_thread_write(fd, '### pseudo pt. time', LOG_HPSI_PSEUDO)
    call timelog_thread_write(fd, '### update time', LOG_HPSI_UPDATE)

    write(fd,*) '### summation time'
    do i=0,NUMBER_THREADS-1
      time = timelog_thread_get(LOG_HPSI_INIT, i)    &
         & + timelog_thread_get(LOG_HPSI_STENCIL, i) &
         & + timelog_thread_get(LOG_HPSI_PSEUDO, i)  &
         & + timelog_thread_get(LOG_HPSI_UPDATE, i)
      write(fd,*) 'tid =',i,':',time,'sec'
    end do

    time = timelog_get(LOG_HPSI_STENCIL)
    write(fd,*) '### global stencil GFLOPS',get_stencil_gflops(time)
    write(fd,*) '### stencil GFLOPS'
    do i=0,NUMBER_THREADS-1
      time    = timelog_thread_get(LOG_HPSI_STENCIL, i)
      cnt     = hpsi_called(i)
      lgflops = get_stencil_gflops(time,cnt)
      write(fd,*) 'tid =',i,': cnt =',cnt,' GFLOPS =',lgflops
      gflops  = gflops + lgflops
    end do
    write(fd,*) '### summation GFLOPS =',gflops
    close(fd)

    print *, ' - summation        :', gflops
  end subroutine
end module opt_variables
